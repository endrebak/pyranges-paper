from __future__ import unicode_literals
from collections import defaultdict
from natsort import natsorted

def supps(w):

    results = _expand(ss[ss.Category == w.category].Function,
                      "{prefix}/benchmark/{function}/{library}/{filetype}/{iteration}_{size}.result",
                      iterations=iterations[:1],
                      sizes=sizes[:1],
                      filetypes="annotation",
                      num_cores=num_cores[:1])

    return results


rule create_test_md:
    input:
        results = supps
    output:
        "supplementary_paper/code_{category}.md"
    run:
        category = wildcards.category
        # result
        code_dict = defaultdict(dict)
        result_dict = defaultdict(dict)
        for f in input.results:

            function, library = f.split("/")[-4:-2]
            library = library.replace("_1", "")

            code_dict[function][library] = category_code[category][library][function]

            result_dict[function][library] = "".join(open(f).readlines())

        outfile = open(output[0], "w+")

        for function, lib_dict in natsorted(code_dict.items()):

            outfile.write("# " + function + "\n\n")

            for library, code in natsorted(lib_dict.items()):
                outfile.write("## " + library + "\n\n")
                outfile.write("#### Code:\n\n")
                outfile.write("```\n" + lib_dict[library] + "\n```\n\n")

                if result_dict[function][library] != "":
                    outfile.write("#### Result:\n\n")
                    outfile.write("```\n" + result_dict[function][library] + "\n```\n\n")

            outfile.write("\\newpage\n")

        outfile.close()


rule create_graph_mds:
    input:
        (f for f in
         expand("{prefix}/benchmark/graphs/{{measure}}_{filetype}_{category}.{extension}", prefix=prefix, filetype=filetypes, category=category_dict, extension=extensions)
         if f.endswith(".png"))
    output:
        "supplementary_paper/{measure}.md"
    run:

        out_handle = open(output[0], "w+")
        out_handle.write("# Timing: PyRanges vs. R GenomicRanges vs. bedtools\n\n")

        graph_dicts = defaultdict(dict)
        for f in input:
            fname = f.split("/")[-1]
            ftype, category = fname.replace(".png", "").split("_")[1:]

            outf = "{}".format(fname)
            command = "cp {} supplementary_paper/{}".format(f, outf)
            print(command)
            shell(command)
            graph_dicts[category][ftype] = fname

        for category, d in graph_dicts.items():

            out_handle.write("### {}\n\n".format(category))

            for ftype in ["annotation", "reads"]:
                out_handle.write("#### {}\n\n".format(ftype))
                img = '<img src="{}" />\n\n'.format(d[ftype])

                out_handle.write(img)

        out_handle.close()


def supps_all(w):

    rows = ss[ss.Function == w.function]
    # print(rows)
    libraries = set(rows.Library)
    # print(libraries)
    libraries = [l.replace("pyranges", "pyranges_1") for l in libraries]

    if w.function == "read_bed" or w.function == "read_bam":
        filetype = "reads"
    else:
        filetype = "annotation"


    results = expand("{prefix}/benchmark/{function}/{library}/{filetype}/0_100000.result",
                     prefix = prefix,
                     filetype = filetype,
                     function = w.function,
                     library=libraries,
                     iterations=iterations[:1],
                     sizes=sizes[:1],
                     filetypes="annotation",
                     num_cores=num_cores[:1])

    return results



def fix_description(desc):

    # from textwrap import wrap
    assert len(desc) < 160, "Description too long (>= 160 chars.)"

    #     constant_desc = """Comparison of the running time and memory usage for PyRanges versus the
    # equivalent libraries in R and/or Python. In the top row the results for GTF
    # data, while the bottom row shows the results for BED data. The left column shows
    # the time usage, while the right column shows the memory usage. Time is measured
    # in log10 seconds, while memory is measured in GigaBytes (GB)."""

    return desc



descriptions = pd.read_csv("supplementaries/descriptions.yaml", sep="\t", header=0)



rule create_per_function_mds:
    input:
        graph = "supplementary_paper/time_memory_together_{function}.png",
        result = supps_all
    output:
        "supplementary_paper/{function}_all.md"
    params:
        description = lambda w: fix_description(descriptions[descriptions.Function == w.function].Description.iloc[0])
    run:
        category = ss[ss.Function == wildcards.function].Category.iloc[0]
        function = wildcards.function

        libraries = ss[ss.Function == function].Library.drop_duplicates().tolist()

        result_handle = open(output[0], "w+")

        result_handle.write("### " + function.capitalize() + "\n\n")

        desc = params.description
        # if wildcards.md_or_pd == "md":
        #     img = '<img src="{}" />\n\n'.format(input.graph.split("/")[-1])
        # else:
        img = '![{desc}]({img} "{fn}")\n\n'.format(img = input.graph.split("/")[-1], fn = function, desc=desc)

        result_handle.write(img)

        result_handle.write("#### Code\n\n")
        for library in libraries:
            result_handle.write("##### " + library + "\n\n")
            result_handle.write("```\n" + category_code[category][library][function] + "\n```\n\n")

        result_handle.write("\pagebreak\n\n")

        # result_handle.write("#### Results\n\n")
        # for library in libraries:
        #     result_handle.write("##### " + library + "\n\n")
        #     f = [f for f in input.result if library.replace("pyranges", "pyranges_1") in f][0]
        #     source = open(f).readlines()
        #     result = "```\n" + "".join(source) + "\n```"
        #     result_handle.write(result + "\n\n")

        result_handle.close()

# "unary": "PyRanges functionality that operates on a single PyRanges. These include functions to sort, cluster, and convert ranges into run-length-encodings (RLE)."

category_descriptions = {"unary": "PyRanges functionality that operates on a single PyRanges object. These include functions to sort, cluster and convert ranges into run length encodings (RLE.)",
                         "binary": "PyRanges functionality that operates on pairs of PyRanges. These functions include functions to find the nearest intervals, find the intersecting intervals, join granges on overlap, set intersect/union and subtract one PyRanges object from another.",
                         "rle": "Arithmetic operations on RLEs. These include add, subtract, divide and multiply.",
                         "tree": "Operations for building and finding overlaps using a tree.",
                         "io": "Functions to read files into PyRanges."}


rule concat_per_category:
    input:
        lambda w: expand("supplementary_paper/{function}_all.md", function=ss[ss.Category == w.category].Function.drop_duplicates())
    output:
        "supplementary_paper/{category}_all.md"
    run:
        oh = open(output[0], "w+")

        content = category_descriptions[wildcards.category]
        oh.write("\pagebreak\n\n## {}\n\n{}\n\n".format(wildcards.category, content))

        for f in input:
            oh.write(open(f).read())

rule concat_categories:
    input:
        expand("supplementary_paper/{category}_all.md", category=category_dict)
    output:
        "supplementary_paper/README.md"
    run:
        oh = open(output[0], "w+")

        title_page = """# PyRanges paper supplementaries

This document shows the time and memory usage for all non-basic functions in
the ecosystem of PyRanges libraries (PyRanges, PyRles and NCLS). They are
compared against their equivalents in Python and R, respectively. The basic
PyRanges functionality is compared against R Bioconductor's GenomicRanges and
pybedtools. The PyRles-functionality is compared against R Bioconductor's
S4Vectors. The NCLS is compared against the intervaltree in the Python
bx-python library. For each function the equivalent code from each library is
shown.\n\n"""
        oh.write(title_page)

        for f in input:
            oh.write(open(f).read())



rule supps_to_pdf:
    input:
        "supplementary_paper/README.md"
    output:
        "supplementary_paper/README.pdf"
    shell:
        "cd supplementary_paper && pandoc --toc --listings -H latex_config.tex -s README.md  -t latex -o README.pdf && cd .."



# rule md_to_pdf:
#     input:
#         "{prefix}/benchmark/supplementaries/{category}.md"
#     output:
#         "{prefix}/benchmark/supplementaries/{category}.pdf"
#     shell:
#         "pandoc -t latex {input[0]} -o {output[0]}"
