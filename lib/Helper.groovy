class Help {

    static def start_info(Map info, String time, String profile) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.1"
        println ""
        if (info.containsKey("fastq")){
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $info.fastq"
        println " Input samples               : $nsamples"
        }
        if (info.containsKey("fasta")){
        println " Input Fasta                 : $info.fasta"
        }
        if (info.containsKey("accessions")){
        println " Input accessions            : $info.accessions"
        }
        println " Reports are found in        : ./reports"
        println " Results are found in        : ./results"
        println " Profile                     : $profile"
        println ""
        println "Starting pipeline at $time"
        println ""

    }

    static void complete_info(nextflow.script.WorkflowMetadata wf) {

        println ""
        println "Pipeline execution summary"
        println "=========================="
        println "Completed at                 : $wf.complete"
        println "Duration                     : $wf.duration"
        println "Success                      : $wf.success"
        println "Work directory               : $wf.workDir"
        println "Exit status                  : $wf.exitStatus"
        println ""

    }

    static def print_help(Map params) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.1"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run pipe.nf"
        println ""
        println "       --fastq                     Path expression to paired-end fastq files. (default: $params.fastq) (default: 'fastq/*_{1,2}.*')"
        println "       "
        println "       Component 'INTEGRITY_COVERAGE_1_1'"
        println "       ----------------------------------"
        println "       --genomeSize_1_1            Genome size estimate for the samples in Mb. It is used to estimate the coverage and other assembly parameters andchecks (default: 1)"
        println "       --minCoverage_1_1           Minimum coverage for a sample to proceed. By default it's setto 0 to allow any coverage (default: 0)"
        println "       "
        println "       Component 'TRIMMOMATIC_1_2'"
        println "       ---------------------------"
        println "       --adapters_1_2              Path to adapters files, if any. (default: 'None')"
        println "       --trimSlidingWindow_1_2     Perform sliding window trimming, cutting once the average quality within the window falls below a threshold (default: '5:20')"
        println "       --trimLeading_1_2           Cut bases off the start of a read, if below a threshold quality (default: 3)"
        println "       --trimTrailing_1_2          Cut bases of the end of a read, if below a threshold quality (default: 3)"
        println "       --trimMinLength_1_2         Drop the read if it is below a specified length  (default: 55)"
        println "       --clearInput_1_2            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'FASTQC_1_3'"
        println "       ----------------------"
        println "       --adapters_1_3              Path to adapters files, if any. (default: 'None')"
        
    }

}

class CollectInitialMetadata {

    public static void print_metadata(nextflow.script.WorkflowMetadata workflow){

        def treeDag = new File("${workflow.projectDir}/.treeDag.json").text
        def forkTree = new File("${workflow.projectDir}/.forkTree.json").text

        def metadataJson = "{'nfMetadata':{'scriptId':'${workflow.scriptId}',\
'scriptName':'${workflow.scriptName}',\
'profile':'${workflow.profile}',\
'container':'${workflow.container}',\
'containerEngine':'${workflow.containerEngine}',\
'commandLine':'${workflow.commandLine}',\
'runName':'${workflow.runName}',\
'sessionId':'${workflow.sessionId}',\
'projectDir':'${workflow.projectDir}',\
'launchDir':'${workflow.launchDir}',\
'startTime':'${workflow.start}',\
'dag':${treeDag},\
'forks':${forkTree}}}"

        def json = metadataJson.replaceAll("'", '"')

        def jsonFile = new File(".metadata.json")
        jsonFile.write json
    }
}