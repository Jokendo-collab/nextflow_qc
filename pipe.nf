#!/usr/bin/env nextflow

import Helper
import CollectInitialMetadata

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

def infoMap = [:]
if (params.containsKey("fastq")){
    infoMap.put("fastq", file(params.fastq).size())
}
if (params.containsKey("fasta")){
    if (file(params.fasta) instanceof LinkedList){
        infoMap.put("fasta", file(params.fasta).size())
    } else {
        infoMap.put("fasta", 1) 
    }
}
if (params.containsKey("accessions")){
    // checks if params.accessions is different from null
    if (params.accessions) {
        BufferedReader reader = new BufferedReader(new FileReader(params.accessions));
        int lines = 0;
        while (reader.readLine() != null) lines++;
        reader.close();
        infoMap.put("accessions", lines)
    }
}

Help.start_info(infoMap, "$workflow.start", "$workflow.profile")
CollectInitialMetadata.print_metadata(workflow)
    

// Placeholder for main input channels
if (params.fastq instanceof Boolean){exit 1, "'fastq' must be a path pattern. Provide value:'$params.fastq'"}
if (!params.fastq){ exit 1, "'fastq' parameter missing"}
IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty { exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

// Placeholder for secondary input channels


// Placeholder for extra input channels


// Placeholder to fork the raw input channel

IN_fastq_raw.set{ trimmomatic_in_1_0 }


IN_genome_size_1_1 = Channel.value(params.genomeSize_1_1)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize__1_1}'")}

IN_min_coverage_1_1 = Channel.value(params.minCoverage_1_1)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage__1_1}'")}

process integrity_coverage_1_1 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_1 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId integrity_coverage_1_1 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(fastq_pair) from trimmomatic_in_1_0
    val gsize from IN_genome_size_1_1
    val cov from IN_min_coverage_1_1
    // This channel is for the custom options of the integrity_coverage.py
    // script. See the script's documentation for more information.
    val opts from Channel.value('')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_encoding'),
        file('*_phred'),
        file('*_coverage'),
        file('*_max_len') into MAIN_integrity_1_1
    file('*_report') optional true into LOG_report_coverage1_1_1
    set sample_id, val("1_1_integrity_coverage"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_integrity_coverage_1_1
set sample_id, val("integrity_coverage_1_1"), val("1_1"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_integrity_coverage_1_1
file ".versions"

    script:
    template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted_1_1 = Channel.create()
MAIN_PreCoverageCheck_1_1 = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity_1_1.choice(LOG_corrupted_1_1, MAIN_PreCoverageCheck_1_1) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
integrity_coverage_out_1_0_dep = Channel.create()
SIDE_phred_1_1 = Channel.create()
SIDE_max_len_1_1 = Channel.create()

MAIN_PreCoverageCheck_1_1
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate(integrity_coverage_out_1_0_dep, SIDE_phred_1_1, SIDE_max_len_1_1){
        a -> [ [a[0], a[1]], [a[0], a[3].text], [a[0], a[5].text]  ]
    }

/** REPORT_COVERAGE - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage_1_1 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_1_1/'

    input:
    file(report) from LOG_report_coverage1_1_1.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}

/** REPORT_CORRUPT - PLUG-IN
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt_1_1 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted_1_1/'

    input:
    val sample_id from LOG_corrupted_1_1.collect{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${sample_id.join(",")} | tr "," "\n" >> corrupted_samples.txt
    """

}


SIDE_phred_1_1.set{ SIDE_phred_1_2 }


// Check sliding window parameter
if ( params.trimSlidingWindow_1_2.toString().split(":").size() != 2 ){
    exit 1, "'trimSlidingWindow_1_2' parameter must contain two values separated by a ':'. Provided value: '${params.trimSlidingWindow_1_2}'"
}
if ( !params.trimLeading_1_2.toString().isNumber() ){
    exit 1, "'trimLeading_1_2' parameter must be a number. Provide value: '${params.trimLeading_1_2}'"
}
if ( !params.trimTrailing_1_2.toString().isNumber() ){
    exit 1, "'trimTrailing_1_2' parameter must be a number. Provide value: '${params.trimTrailing_1_2}'"
}
if ( !params.trimMinLength_1_2.toString().isNumber() ){
    exit 1, "'trimMinLength_1_2' parameter must be a number. Provide value: '${params.trimMinLength_1_2}'"
}

IN_trimmomatic_opts_1_2 = Channel.value([params.trimSlidingWindow_1_2,params.trimLeading_1_2,params.trimTrailing_1_2,params.trimMinLength_1_2])
IN_adapters_1_2 = Channel.value(params.adapters_1_2)

clear = params.clearInput_1_2 ? "true" : "false"
checkpointClear_1_2 = Channel.value(clear)

process trimmomatic_1_2 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_2 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId trimmomatic_1_2 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    publishDir "results/trimmomatic_1_2", pattern: "*.gz"

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair), phred from integrity_coverage_out_1_0_dep.join(SIDE_phred_1_2)
    val trim_range from Channel.value("None")
    val opts from IN_trimmomatic_opts_1_2
    val ad from IN_adapters_1_2
    val clear from checkpointClear_1_2

    output:
    set sample_id, "${sample_id}_*trim.fastq.gz" into trimmomatic_out_1_0
    file 'trimmomatic_report.csv'
    set sample_id, val("1_2_trimmomatic"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_trimmomatic_1_2
set sample_id, val("trimmomatic_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_trimmomatic_1_2
file ".versions"

    script:
    template "trimmomatic.py"

}



IN_adapters_1_3 = Channel.value(params.adapters_1_3)

process fastqc2_1_3 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_1_3 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir "reports/fastqc_1_3/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from trimmomatic_out_1_0
    val ad from IN_adapters_1_3

    output:
    set sample_id, file(fastq_pair), file('pair_1*'), file('pair_2*') into MAIN_fastqc_out_1_3
    file "*html"
    set sample_id, val("1_3_fastqc2"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc2_1_3
set sample_id, val("fastqc2_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc2_1_3
file ".versions"

    script:
    template "fastqc.py"
}


process fastqc2_report_1_3 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_1_3 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc_1_3/run_2/', pattern: '*summary.txt', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_1_3
    val opts from Channel.value("")

    output:
    set sample_id, file(fastq_pair), '.status' into MAIN_fastqc_report_1_3
    file "*_status_report" into LOG_fastqc_report_1_3
    file "${sample_id}_*_summary.txt" optional true
    set sample_id, val("1_3_fastqc2_report"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc2_report_1_3
set sample_id, val("fastqc2_report_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc2_report_1_3
file ".versions"

    script:
    template "fastqc_report.py"

}


process compile_fastqc_status2_1_3 {

    publishDir 'reports/fastqc_1_3/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_1_3.collect()

    output:
    file 'FastQC_2run_report.csv'

    """
    echo Sample, Failed? >> FastQC_2run_report.csv
    cat $rep >> FastQC_2run_report.csv
    """

}

fastqc_out_1_1 = Channel.create()

MAIN_fastqc_report_1_3
        .filter{ it[2].text == "pass" }
        .map{ [it[0], it[1]] }
        .into(fastqc_out_1_1)




/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { sample_id }
    publishDir "pipeline_status/$task_name"

    input:
    set sample_id, task_name, status, warning, fail, file(log) from STATUS_integrity_coverage_1_1.mix(STATUS_trimmomatic_1_2,STATUS_fastqc2_1_3,STATUS_fastqc2_report_1_3)

    output:
    file '*.status' into master_status
    file '*.warning' into master_warning
    file '*.fail' into master_fail
    file '*.log'

    """
    echo $sample_id, $task_name, \$(cat $status) > ${sample_id}_${task_name}.status
    echo $sample_id, $task_name, \$(cat $warning) > ${sample_id}_${task_name}.warning
    echo $sample_id, $task_name, \$(cat $fail) > ${sample_id}_${task_name}.fail
    echo "\$(cat .command.log)" > ${sample_id}_${task_name}.log
    """
}

process compile_status_buffer {

    input:
    file status from master_status.buffer( size: 5000, remainder: true)
    file warning from master_warning.buffer( size: 5000, remainder: true)
    file fail from master_fail.buffer( size: 5000, remainder: true)

    output:
    file 'master_status_*.csv' into compile_status_buffer
    file 'master_warning_*.csv' into compile_warning_buffer
    file 'master_fail_*.csv' into compile_fail_buffer

    """
    cat $status >> master_status_${task.index}.csv
    cat $warning >> master_warning_${task.index}.csv
    cat $fail >> master_fail_${task.index}.csv
    """
}

process compile_status {

    publishDir 'reports/status'

    input:
    file status from compile_status_buffer.collect()
    file warning from compile_warning_buffer.collect()
    file fail from compile_fail_buffer.collect()

    output:
    file "*.csv"

    """
    cat $status >> master_status.csv
    cat $warning >> master_warning.csv
    cat $fail >> master_fail.csv
    """

}


/** Reports
Compiles the reports from every process
*/
process report {

    tag { sample_id }

    input:
    set sample_id,
            task_name,
            pid,
            report_json,
            version_json,
            trace from REPORT_integrity_coverage_1_1.mix(REPORT_trimmomatic_1_2,REPORT_fastqc2_1_3,REPORT_fastqc2_report_1_3)

    output:
    file "*" optional true into master_report

    """
    prepare_reports.py $report_json $version_json $trace $sample_id $task_name 1 $pid $workflow.scriptId $workflow.runName
    """

}


process compile_reports {

    publishDir "pipeline_report/", mode: "copy"

    if ( params.reportHTTP != null ){
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH;"
        afterScript "metadata_POST.sh $params.projectId $params.pipelineId 0 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId 0 \"$params.platformSpecies\""
    }

    input:
    file report from master_report.collect()
    file forks from Channel.fromPath("${workflow.projectDir}/.forkTree.json")
    file dag from Channel.fromPath("${workflow.projectDir}/.treeDag.json")
    file js from Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")

    output:
    file "pipeline_report.json"
    file "pipeline_report.html"
    file "src/main.js"

    script:
    template "compile_reports.py"
}

workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
