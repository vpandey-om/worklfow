process RENDER_REPORT {
    tag "report"
    label 'process_low'

    input:
    path inputs

    output:
    path "survom_report.html", emit: html

    script:
    """
    echo "<html><body><h1>SurvOm report</h1></body></html>" > survom_report.html
    """
}
