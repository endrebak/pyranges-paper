from collections import defaultdict
from natsort import natsorted

def supps(w):

    # print("w: ", w)
    results = _expand(ss[ss.Category == w.category].Function,
                      "{prefix}/benchmark/{function}/{library}/{filetype}/{iteration}_{size}.result",
                      iterations=iterations[:1],
                      sizes=sizes[:1],
                      filetypes="annotation",
                      num_cores=num_cores[:1])

    print(results)
    return results


rule create_test_md:
    input:
        results = supps
    output:
        "{prefix}/benchmark/supplementaries/{category}.md"
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


rule md_to_pdf:
    input:
        "{prefix}/benchmark/supplementaries/{category}.md"
    output:
        "{prefix}/benchmark/supplementaries/{category}.pdf"
    shell:
        "pandoc -t latex {input[0]} -o {output[0]}"
