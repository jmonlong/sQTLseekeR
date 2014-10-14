makeClusterFunctionsAdaptive <- function (template.file){
    library(stringr)
    template = cfReadBrewTemplate(template.file)
    submitJob = function(conf, reg, job.name, rscript, log.file,
        job.dir, resources, arrayjobs) {
        outfile = cfBrewTemplate(conf, template, rscript, "pbs")
        res =  BatchJobs:::runOSCommandLinux("qsub", outfile, stop.on.exit.code = FALSE)
        max.jobs.msg = "Maximum number of jobs already in queue"
        output = collapse(res$output, sep = "\n")
        if (grepl(max.jobs.msg, output, fixed = TRUE)) {
            makeSubmitJobResult(status = 1L, batch.job.id = NA_character_,
                                msg = max.jobs.msg)
        }
        else if (res$exit.code > 0L) {
            cfHandleUnknownSubmitError("qsub", res$exit.code,
                                       res$output)
        }
        else {
            makeSubmitJobResult(status = 0L, batch.job.id = str_extract(output,"\\d+"))
        }
    }
    killJob = function(conf, reg, batch.job.id) {
        cfKillBatchJob("qdel", batch.job.id)
    }
    listJobs = function(conf, reg) {
        BatchJobs:::runOSCommandLinux("qstat", "-u $USER")$output
    }
    getArrayEnvirName = function() "TASKID"
    makeClusterFunctions(name = "Adaptive", submitJob = submitJob,
                         killJob = killJob, listJobs = listJobs, getArrayEnvirName = getArrayEnvirName)
}
