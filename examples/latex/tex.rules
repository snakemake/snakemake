ruleorder:  tex2pdf_with_bib > tex2pdf_without_bib

rule tex2pdf_with_bib:
    input:
        '{name}.tex',
        '{name}.bib'
    output:
        '{name}.pdf'
    shell:
        """
        pdflatex {wildcards.name}
        bibtex {wildcards.name}
        pdflatex {wildcards.name}
        pdflatex {wildcards.name}
        """

rule tex2pdf_without_bib:
    input:
        '{name}.tex'
    output:
        '{name}.pdf'
    shell:
        """
        pdflatex {wildcards.name}
        pdflatex {wildcards.name}
        """

rule texclean:
    shell:
        "rm -f  *.log *.aux *.bbl *.blg *.synctex.gz"
