#!/usr/bin/php
<?php

// Variables used in the try/catch block.
$result = null;
$jobStatus = null;
$stdError = null;
$fileSummaries;
$totals;

// Variables used by the command below.
$proposalsPath = "/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load/testData/msl39v3/proposals_msl39v3_crash1";
$resultsPath = "/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load/testResultsDocker/msl39v3/proposals_msl39v3_crash1";
$scriptName = "curtish/ictv_proposal_processor:v2.1.096009a"; // crash
$scriptName = "curtish/ictv_proposal_processor:v2.2.acb51d8"; // no crash
echo("proposalsPath: ".$proposalsPath."\n");
echo("resultsPath:   ".$resultsPath."\n");
echo("scriptName:    ".$scriptName."\n");

// The directory where the command will be run.
$workingDirectory = "./";

// An indexed array where the key represents the descriptor number and the value represents how 
// PHP will pass that descriptor to the child process. 0 is stdin, 1 is stdout, while 2 is stderr.
$descriptorSpec = array(
   0 => array("pipe", "r"), // Read from stdin (not used)
   1 => array("pipe", "w"), // Write to stdout
   2 => array("pipe", "w")  // Write to stderr
);
        
// Generate the command to be run.
$command = "docker run ".
   "-v \"{$proposalsPath}:/proposalsTest\":ro ".
   "-v \"{$resultsPath}:/results\" ".
   $scriptName." ".
   "/merge_proposal_zips.R -v ";

echo("Command: ".$command."\n");

try {
   // See https://www.php.net/manual/en/function.proc-open.php for details.
   $process = proc_open($command, $descriptorSpec, $pipes, $workingDirectory);
   
   echo("Process opened: ".$process."\n");

   if (is_resource($process)) {

	 echo("Process is_resource: TRUE\n");

         // $pipes now looks like this:
         // 0 => writeable handle connected to child stdin
         // 1 => readable handle connected to child stdout
         // 2 => writeable handle connected to child stderr

         // Note: We're not using the stdin pipe.

         // Get stdout
         $result = stream_get_contents($pipes[1]);
	 echo("stdout: \n".$result."\n/stdout\n");
         fclose($pipes[1]);

         // Get stderror
         $stdError = stream_get_contents($pipes[2]);
	 echo("stderr: \n".$stdError."\n/stderr\n");
         fclose($pipes[2]);

         // It is important that you close any pipes before calling proc_close in order to avoid a deadlock.
         $exitCode = proc_close($process);
	 echo("exitCode: ".$exitCode."\n");

         // In this case "pending complete validation".
         $jobStatus = "pending";
	 echo("jobStatus: ".$jobStatus."\n");

   } else {
         $jobStatus = "crashed";
	 echo("jobStatus: ".$jobStatus."\n");
         $stdError = "Process is not a resource";
   }

   if ($jobStatus != "crashed") {

      echo("not crashed: jobStatus: ".$jobStatus."\n");
      // Curtis: The commented code below is what actually happens in the proposal service.

      // Parse the summary TSV file for proposal filenames and their status counts (file summaries).
      $fileSummaries = ProposalFileSummary::getFileSummaries($resultsPath);
      echo("fileSummaries: \n".$fileSummaries."\n/fileSummaries\n");

      // If no file summaries were found, return a job status of "crashed".
      if (!$fileSummaries || sizeof($fileSummaries) < 1) { 
	      $jobStatus = JobStatus::$crashed; 
	      echo("no fileSummaries found: set jobStatus=".$jobStatus."\n");
      }
   }
} 
catch (Exception $e) {
   echo("Exception: ".$e->getMessage()."\n");

   $jobStatus = "crashed";

   if ($e) { 
      if (isset($stdError) && $stdError !== '') { $stdError = $stdError . "; "; }
      $stdError = $stdError.$e->getMessage(); 
   }

   echo("The following error occurred: ".$stdError);
}
