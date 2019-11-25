<?php

require_once('AlignmentClass.php');

//For each word match, extend alignment in both directions to find
//alignments that score greater than score threshold S.
function global_alignment($query_seq_upper, $protein_subject, $matching_words_indices, $score_threshold, $word_size, $matrix, $gap_costs) {
    $query_match_index = 0;
    $db_match_index = 0;
    $index = 0;
    $optimal_alignments = array();
    $sequence_alignments = array();
    $unique_alignments = array();
    $matching_words_size = count($matching_words_indices);
    while ($index < $matching_words_size) {
        $query_match_size = count($matching_words_indices[$index][0]);
        $db_match_size = count($matching_words_indices[$index][1]);
        for ($query_index = 0; $query_index < $query_match_size; $query_index++) {
            for ($db_index = 0; $db_index < $db_match_size; $db_index++) {
                $query_match_index = $matching_words_indices[$index][0][$query_index];
                $db_match_index = $matching_words_indices[$index][1][$db_index];
                // get optimal alignments by needleman algorithm global alignments
                // global alignments over the whole sequence
                $alignments = new AlignmentClass($query_seq_upper, $protein_subject, $query_match_index, $db_match_index, $score_threshold, $word_size, $gap_costs, $matrix);
                // for each alignment, add to the collection sequence alignments
                $sequence_alignments = $alignments->get_global_alignments();
                // push $sequence_alignments into $optimal_alignments
                array_push($optimal_alignments, $sequence_alignments);
            }
        }
        $index++;
    } // end while loop
    // sort the aligments by most significant on top
    $optimal_alignments = mostSignificantAlignments($optimal_alignments);

    // get unique alignments only
    $unique_alignments = getUnique_alignments($optimal_alignments);

    return $unique_alignments;
}

// get unique alignments only
function getUnique_alignments($alignments) {
    $key_array = array();
    $temp_array = array();
    $alignment_size = count($alignments);
    for ($i = 0; $i < $alignment_size; $i++) {
        //if have not seen alignment in array then insert it.
        if (!in_array($alignments[$i]['optimal']['match_score'], $key_array, true)) {
            array_push($temp_array, $alignments[$i]);
            array_push($key_array, $alignments[$i]['optimal']['match_score']);
        }
    }

    return $temp_array;
}

// sort the aligments by most significant on top
function mostSignificantAlignments($alignments) {
    $current = 0;
    $alignments_size = count($alignments);
    while ($current < $alignments_size) {
        // swap the most significant aligments index to the top
        for ($next = $current + 1; $next < $alignments_size; $next++) {
            $next_aligment_index = $alignments[$next]['optimal']['match_score'];
            $current_aligment_index = $alignments[$current]['optimal']['match_score'];
            // if next alignment score is higher than current alignment
            // swap their position.
            if ($next_aligment_index > $current_aligment_index) {
                $most_significant = $alignments[$next];
                $alignments[$next] = $alignments[$current];
                $alignments[$current] = $most_significant;
            }
        }
        $current++;
    }
    return $alignments;
}

// score threshold counter
function scoreThreshold($sub_query, $sub_subject, $thr) {
    $compare_sequence = strcmp($sub_query, $sub_subject);
    // if not equal then increase threshold counter
    if ($compare_sequence != 0) {
        $thr++;
    }
    // else decrement counter if not zero
    elseif ($thr != 0) {
        $thr--;
    }

    return $thr;
}

// find locations of matching words in database sequences
// Compare the query word list to the database and identify exact matches.
// save the indices
function matching_words($query_words_list, $protein_subject_words_list) {
    $index = 0;
    $query_size = count($query_words_list);
    $matching_words_indices = array();

    // go through the whole query words list
    while ($index < $query_size) {
        // key — Fetch a key from an array
        $query_word_key = key($query_words_list[$index]);
        $db_index = 0;
        $db_size = count($protein_subject_words_list);

        // go through each database word list
        while ($db_index < $db_size) {
            $word_index_match = array();
            $db_word_key = key($protein_subject_words_list[$db_index]);
            // Returns < 0 if str1 is less than str2; > 0 if str1 is greater 
            // than str2, and 0 if they are equal.
            $word_match = strcmp($query_word_key, $db_word_key);
            // if query word match with database then push it to array
            if ($word_match == 0) {
                $query_word_index = $query_words_list[$index][$query_word_key];
                $db_word_index = $protein_subject_words_list[$db_index][$db_word_key];
                $word_index_match = array($query_word_index, $db_word_index);
                array_push($matching_words_indices, $word_index_match);
            }
            $db_index++;
        } // end db
        $index++;
    }

    return $matching_words_indices;
}

// breaking sequence into words and make a word list out of it.
// key : word
// value : indices (location)
function words_list($sequence, $word_size) {
    $word_list = array();
    $sequence_size = strlen($sequence);
    $index = 0;

    //go through the sequence to get each word
    while ($index < $sequence_size) {
        $flag = FALSE;
        $current_word = $index + $word_size - 1;
        $word_list_count = count($word_list);
        //substr ( string $string , int $start [, int $length ] ) : string
        //Returns the portion of string specified by the start and length parameters.
        $new_word = substr($sequence, $index, $word_size);

        // base case
        if ($current_word >= $sequence_size) {
            break;
        }

        $list_index = 0;
        while ($list_index < $word_list_count) {
            // check if word_list_index isset then push it into the array
            // already seen word
            if (isset($word_list[$list_index][$new_word])) {
                array_push($word_list[$list_index][$new_word], $index);
                $flag = TRUE;
                break;
            }
            $list_index++;
        }

        // not seen word in word list, add word
        if (!$flag) {
            $never_seen = array($new_word => array($index));
            array_push($word_list, $never_seen);
            $flag = FALSE;
        }
        $index++;
    }

    return $word_list;
}

?>