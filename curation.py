#!/usr/bin/python
import argparse
import sys
import suds
import re
from suds.client import Client
import numpy as np
import pandas as pd
import time
import requests
from urllib.request import urlretrieve
from urllib.parse import quote
import socket
from collections import Counter
import plotly.express as px

start_time = time.time()
server = "https://rest.ensembl.org"
ext = "/vep/human/hgvs/"
URL = 'https://mutalyzer.nl/services/?wsdl'
c = Client(URL, cache=None)
o = c.service


# ---parse commandline arguments---
def parse_args():
    parser = argparse.ArgumentParser(description="Script for curation of BRCA database")
    parser.add_argument("-i", "--input", help="path to input file")
    parser.add_argument("-o", "--output", help="path to new output file, default: overwrite old file")
    parser.add_argument("-w", "--warnings", help="display warnings and errors, default=false", default=False)
    parser.add_argument("-s", "--stat", help="display statistic about curation, default=false", default=False)
    return vars(parser.parse_args())


def figures_to_html(figs, filename="stat.html"):
    dashboard = open(filename, 'w')
    dashboard.write("<html><head></head><body>" + "\n")
    for fig in figs:
        inner_html = fig.to_html().split('<body>')[1].split('</body>')[0]
        dashboard.write(inner_html)
    dashboard.write("</body></html>" + "\n")


# -----displays statistic over errors that occured during curation and they're frequency----#
def statistic(errors_genomic, errors_protein, errors_conseq):
    eg_count = errors_genomic["error"].value_counts().rename_axis('error').reset_index(name='counts')
    ep_count = errors_protein["error"].value_counts().rename_axis('error').reset_index(name='counts')
    ec_count = errors_conseq["error"].value_counts().rename_axis('error').reset_index(name='counts')
    tit = "unable to complete genomic alterations: " + str(len(errors_genomic.index)) + " from " + str(
         errors_genomic["total"][0])
    print(tit)
    fig1 = px.pie(eg_count, values='counts', names='error', title=tit)
    tit = "unable to complete protein alterations: " + str(len(errors_protein.index)) + " from " + str(
         errors_protein["total"][0])
    fig2 = px.pie(ep_count, values='counts', names='error', title=tit)
    tit = "unable to complete protein alterations: " + str(len(errors_conseq.index)) + " from " + str(
         errors_conseq["total"][0])
    fig3 = px.pie(ec_count, values='counts', names='error', title=tit)
    same = (set(errors_conseq["index"]).intersection(errors_protein["index"])).intersection(errors_genomic["index"])
    x = [errors_genomic["error"][list(errors_genomic["index"]).index(s)] for s in same]
    y = [errors_protein["error"][list(errors_protein["index"]).index(s)] for s in same]
    z = [errors_conseq["error"][list(errors_conseq["index"]).index(s)] for s in same]
    same_err = [e for i, e in enumerate(x) if (e in y[i]) & (e in z[i])]
    same_err = pd.DataFrame({"error": list(set(same_err)), "counts": [same_err.count(x) for x in set(same_err)]})
    fig4 = px.pie(same_err, values='counts', names='error',
                  title='unable to complete GA, PA and consequence in same row with same error code')
    figures_to_html([fig1, fig2, fig3, fig4])


# ----add missing genomic alterations based on cDNA coordinates via Mutalyzer numberConversion-----
def genomicAlt(dat, warning, stat, ref_seqs):
    indices = []
    errors = []
    missing_dna = dat.loc[dat['GENOMISCHE_ALTERATION'].isnull()]
    for row in missing_dna.itertuples():
        if row.GEN not in ref_seqs:
            if warning:
                print("unkown Gene")
            if stat:
                errors.append("unkown Gene")
                indices.append(row.Index)
            continue
        variant = ref_seqs[row.GEN] + ":" + row.DNA_HGVS
        r = o.checkSyntax(variant)
        if warning:
            if r.messages:
                for m in r.messages.SoapMessage:
                    print('Message (%s): %s' % (m.errorcode, m.message))
        if r.valid:
            genom = o.numberConversion("hg19", variant)
            if genom is None or not genom[0]:
                if warning:
                    print("no conversion possible of variant: ", variant)
                if stat:
                    errors.append("unable to convert")
                    indices.append(row.Index)
            else:
                genom_variant = re.sub(r'.*g', 'g', genom[0][0])
                dat.loc[row.Index, 'GENOMISCHE_ALTERATION'] = genom_variant
        else:
            if warning:
                print('variant not valid!', " ", variant)
            if stat:
                if row.DNA_HGVS[0:2] != "c.":
                    errors.append("missing cDNA")
                else:
                    errors.append("unable to parse cDNA, syntax error")
                indices.append(row.Index)
    errors_genomic = pd.DataFrame(
        {"index": indices, "error": errors, "total": ([len(missing_dna.index)] * len(errors))})
    return dat, errors_genomic


# -----add missing protein alterations using genomic coordinates or cDNA and Mutalyzer----#
def proteinAlt(dat, warning, stat, ref_seqs, acc_nums):
    indices = []
    errors = []
    missing_prot = dat.loc[(dat['PROTEIN_HGVS'].isnull()) | (dat['PROTEIN_HGVS'] == "<keine Angabe>") | (
            dat['PROTEIN_HGVS'] == "unbekannt")]
    ecode = {"EPARSE": "unable to parse cDNA, syntax error", "EREF": "found different base at refseq",
             'ENOINTRON': "Intronic position given for non-genomic reference sequence.",
             'ENOTIMPLEMENTED': 'Only insertion of sequence/range is implemented.',
             'EARGLEN': 'length of variant differed from range.', "ERANGE": "end pos is smaller that begin pos"}
    for row in missing_prot.itertuples():
        if row.GEN not in ref_seqs:
            if warning:
                print("unkown Gene")
            if stat:
                errors.append("unkown Gene")
                indices.append(row.Index)
            continue
        variant = acc_nums[row.GEN] + ":" + str(row.GENOMISCHE_ALTERATION)
        if str(row.GENOMISCHE_ALTERATION) == "nan":
            variant = ref_seqs[row.GEN] + ":" + str(row.DNA_HGVS)
        try:
            r = o.runMutalyzerLight(variant)
        except socket.timeout:
            if warning:
                print("timeout error")
            if stat:
                errors.append("timeout error")
                indices.append(row.Index)
            continue

        prot = r["proteinDescriptions"]
        if prot is None or not prot[0]:
            if warning:
                print("no protein description found for ", variant)
            if stat:
                if (str(row.DNA_HGVS)[0:2] != "c.") & (str(row.GENOMISCHE_ALTERATION)[0:2] != "g."):
                    errors.append("missing cDNA/DNA")
                else:
                    errors.append(",".join(set(
                        [ecode[m.errorcode] if m.errorcode in ecode else m.message for m in r.messages.SoapMessage])))
                indices.append(row.Index)
        else:
            var = [re.sub(r'.*p', 'p', s) for s in prot[0] if row.GEN in s]
            if len(set(var)) == 1:
                protein = var[0]
                dat.loc[row.Index, 'PROTEIN_HGVS'] = protein

            else:
                if warning:
                    print("ambiguous protein description")
                if stat:
                    errors.append("ambiguous protein description")
                    indices.append(row.Index)
    errors_protein = pd.DataFrame(
        {"index": indices, "error": errors, "total": ([len(missing_prot.index)] * len(errors))})
    return dat, errors_protein


# ------ add missing consequences using Ensemble VEP -------
def consequence(dat, warning, stat):
    indices = []
    errors = []
    missing_conseq = dat.loc[
        (dat['KONSEQUENZ'].isnull()) | (dat['KONSEQUENZ'] == "keine Angabe") | (dat['KONSEQUENZ'] == "unbekannt") | (
                dat["KONSEQUENZ"] == "unklar")]
    for row in missing_conseq.itertuples():
        variant = ext + row.GEN + ":" + row.DNA_HGVS
        u = server + variant
        try:
            r = requests.get(u, headers={"Content-Type": "application/json"})
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            if warning:
                print(r.text)
            if stat:
                indices.append(row.Index)
                if str(row.DNA_HGVS)[0:2] != "c.":
                    errors.append("missing cDNA/DNA")
                elif ("can't interpret" in r.text) | ("parse" in r.text):
                    errors.append("unable to parse cDNA, syntax error")
                    continue
                elif "does not match reference allele" in r.text:
                    errors.append("found different base at refseq")
                    continue
                else:
                    errors.append(r.text)
        if r.ok:
            decoded = r.json()
            conseq = ",".join(decoded[0]["transcript_consequences"][0]["consequence_terms"])
            dat.loc[row.Index, 'KONSEQUENZ'] = conseq
    errors_conseq = pd.DataFrame(
        {"index": indices, "error": errors, "total": ([len(missing_conseq.index)] * len(errors))})
    return dat, errors_conseq


def main():
    args = parse_args()
    filename = args["input"]
    warnings = args["warnings"]
    stat = args["stat"]
    dat = pd.read_excel(filename, 0)
    print("open")
    genes = dat["GEN"].unique()
    ref_seq = []
    acc_num = []
    for g in genes:
        try:
            o.getTranscriptsByGeneName("hg19", g)
            o.getGeneLocation(g, "hg19")
        except Exception as e:
            if warnings:
                print(e)
            continue
        ref_seq.append(o.getTranscriptsByGeneName("hg19", g)[0][0])
        acc_num.append(o.getGeneLocation(g, "hg19")["chromosome_accession"])
    # ref_seq = [o.getTranscriptsByGeneName("hg19", g)[0][0] for g in genes]
    # acc_num = [o.getGeneLocation(g, "hg19")["chromosome_accession"] for g in genes]
    acc_nums = dict(zip(genes, acc_num))
    ref_seqs = dict(zip(genes, ref_seq))

    dat, eg = genomicAlt(dat, warnings, stat, ref_seqs)
    dat, ep = proteinAlt(dat, warnings, stat, ref_seqs, acc_nums)
    dat, ec = consequence(dat, warnings, stat)

    if "output" in args:
        filename = args["output"]
    dat.to_excel(filename)

    if stat:
        statistic(eg, ep, ec)

    if warnings:
        # ------------ print time ------------
        seconds = time.time() - start_time
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        print(h, ":", m, ":", s)


if __name__ == "__main__":
    main()
