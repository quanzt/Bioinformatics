<?php
$doc_root = filter_input(INPUT_SERVER, 'DOCUMENT_ROOT');

// Improved way to set the include path to the project root
// Works even if the project is redeployed at another
// level in the web server's filesystem
$dirs = explode(DIRECTORY_SEPARATOR, __DIR__);
array_pop($dirs); // remove last element
$project_root = implode('/',$dirs) . '/';
ini_set('max_execution_time', 1000); 
ini_set('display_errors', '0');
ini_set('log_errors', 1);
// the following file needs to exist, be accessible to apache
// and writable (on Linux: chmod 777 php-errors.log,
// Windows defaults to writable)
// Use an absolute file path to create just one log for the web app
ini_set('error_log', $project_root . 'php-errors.log');

require('./models/global_alignment_db.php');

// landing page: welcome page
$action = filter_input(INPUT_POST, 'action');
if ($action == NULL) {
    $action = filter_input(INPUT_GET, 'action');
    if ($action == NULL) {
        $action = 'blastp_welcome';
    }
}

if ($action == 'blastp_welcome') {
    include('blastp_welcome.php');
}

// get alignment result
elseif ($action == 'getBlastp_result') {
    try {
        // assign input into data field
        $query_seq = filter_input(INPUT_POST, 'query_seq');
        $query_seq_upper = strtoupper($query_seq);
        $protein = filter_input(INPUT_POST, 'protein');
        $score_threshold = filter_input(INPUT_POST, 'score_threshold', FILTER_VALIDATE_INT);
        $word_size = filter_input(INPUT_POST, 'word_size', FILTER_VALIDATE_INT);
        $matrix = filter_input(INPUT_POST, 'matrix');
        $gap_costs = filter_input(INPUT_POST, 'gap_costs', FILTER_VALIDATE_INT);
    } catch (Exception $alignment) {
        include('../errors/error.php');
        exit();
    }

    //call functions to get alignments result
    try {
        // get the protein subject sequence file to use it as database sequence
        $protein_subj_file = file($protein, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
        $protein_subject = $protein_subj_file[0];
        // Break query into words
        $query_words_list = words_list($query_seq_upper, $word_size);
        // Break database sequence into words
        $protein_subject_words_list = words_list($protein_subject, $word_size);
        // find locations of matching words in database sequences
        $matching_words_indices = matching_words($query_words_list, $protein_subject_words_list);
        // get alignments by global alignment algorithm.
        $alignments = global_alignment($query_seq_upper, $protein_subject, $matching_words_indices, $score_threshold, $word_size, $matrix, $gap_costs);
    } catch (Exception $alignment) {
        include('../errors/error.php');
        exit();
    }
    include 'alignments_view.php';
}

?>