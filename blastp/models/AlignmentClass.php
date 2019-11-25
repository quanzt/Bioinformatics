<?php

require_once('NeedlemanClass.php');

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

class AlignmentClass {

    private $query_sequence = '';
    private $protein_subject = '';
    private $query_match_index = 0;
    private $db_match_index = 0;
    private $score_threshold = 0;
    private $word_size = 0;
    private $gap_costs = 0;
    private $matrix = array();
    private $current = 0;
    private $statistical_threshold = 0;
    private $query_seq_extend_right = '';
    private $db_seq_extend_right = '';
    private $query_seq_extend_left = '';
    private $db_seq_extend_left = '';
    private $terminate_extension_right = FALSE;
    private $terminate_left_extension = FALSE;
    private $query_seed = '';
    private $db_seed = '';
    private $score_array = array();
    private $optimal_match_score = 0;
    private $query_score_extend_left = '';
    private $query_score_extend_right = '';
    private $db_score_extend_left = '';
    private $db_score_extend_right = '';
    private $query_needleman_alignment = '';
    private $db_needleman_alignment = '';
    private $match_needleman_alignment = '';

    public function __construct($query_seq_upper, $protein_subject, $query_match_index, $db_match_index, $score_threshold, $word_size, $gap_costs, $matrix) {

        $this->query_sequence = $query_seq_upper;
        $this->protein_subject = $protein_subject;
        $this->query_match_index = $query_match_index;
        $this->db_match_index = $db_match_index;
        $this->score_threshold = $score_threshold;
        $this->word_size = $word_size;
        $this->gap_costs = $gap_costs;
        $this->matrix = $matrix;
    }

    public function get_global_alignments() {
        $needleman_alignments = array();
        // start at the middle of the sequence
        //substr ( string $string , int $start [, int $length ] ) : string
        //Returns the portion of string specified by the start and length parameters.
        $word_median = substr($this->query_sequence, $this->query_match_index, $this->word_size);
        $this->query_seed = $word_median;
        $this->db_seed = $word_median;
        // needleman object
        $needlemanAlgorithm = new NeedlemanWunsch($this->gap_costs, $this->matrix);
        // Stop extension if alignment score fall below score threshold
        while ($this->statistical_threshold < $this->score_threshold) {
            
            // extend from query seed and the db seed to right by one letter/base
            $query_next_start_right = $this->query_match_index + $this->word_size + $this->current;
            $db_next_start_right = $this->db_match_index + $this->word_size + $this->current;
            $query_subproblem = substr($this->query_sequence, $query_next_start_right, 1);
            $db_subproblem = substr($this->protein_subject, $db_next_start_right, 1);
            $flag = FALSE;
            // if two word matches as anchors to build an optimal alignments from the subproblem
            if (!empty($query_subproblem)) {
                if (!empty($db_subproblem)) {                    
                    // use needlman algorithm to get optimal sequence
                    $needleman_alignments = $this->get_needleman($needlemanAlgorithm, $query_subproblem, $db_subproblem, $this->query_seed, $this->db_seed, $flag);
                    $this->set_optimal($needleman_alignments);
                } else {
                    $this->terminate_extension_right = TRUE;
                }
            } else {
                $this->terminate_extension_right = TRUE;
            }

            
            // extend from query seed and the db seed to left by one letter/base
            $query_next_start_left = $this->query_match_index - $this->current - 1;
            $db_next_start_left = $this->db_match_index - $this->current - 1;
            $query_subproblem = substr($this->query_sequence, $query_next_start_left, 1);
            $db_subproblem = substr($this->protein_subject, $db_next_start_left, 1);
            // if two word matches as anchors to build an optimal alignments from the subproblem
            if ($query_next_start_left >= 0) {
                if ($db_next_start_left >= 0) {
                    $flag = TRUE;
                    // use needlman algorithm to get optimal sequence
                    $needleman_alignments = $this->get_needleman($needlemanAlgorithm, $query_subproblem, $db_subproblem, $this->query_seed, $this->db_seed, $flag);
                    $this->set_optimal($needleman_alignments);
                } else {
                    $this->terminate_left_extension = TRUE;
                }
            } else {
                $this->terminate_left_extension = TRUE;
            }

            // check if extensions stopped
            if ($this->terminate_extension_right && $this->terminate_left_extension) {
                break;
            } // break out while loop

            $this->current++;
        } // end while loop
        //
        // set and return alignments to global_alignment_db.php which return to view
        $sequence_alignment = $this->setSequence_alignment();

        return $sequence_alignment;
    }

    // set the sequence alignment and return to view
    private function setSequence_alignment() {
        $sequence_alignment = array(
            'query_start_index' => strpos($this->query_sequence, str_replace('-', '', $this->query_needleman_alignment)),
            'query_end_index' => $this->query_match_index + $this->word_size + strlen($this->query_score_extend_right) - 1,
            'db_start_index' => strpos($this->protein_subject, str_replace('-', '', $this->db_needleman_alignment)),
            'db_end_index' => $this->db_match_index + $this->word_size + strlen($this->db_score_extend_right) - 1,
            'optimal' => $this->match_needleman_alignment);
        return $sequence_alignment;
    }

    // use needleman algorithm to get optimal alignment
    private function get_needleman($needlemanAlgorithm, $query_subproblem, $db_subproblem, $query_seed, $db_seed, $flag) {

        // if query and db subproblem not equal 
        // then increase threshold counter else decrement counter
        $this->statistical_threshold = scoreThreshold($query_subproblem, $db_subproblem, $this->statistical_threshold);

        // extension right sequence
        if ($flag == FALSE) {
            $this->query_seq_extend_right .= $query_subproblem;
            $query_seed .= $query_subproblem;
            $this->db_seq_extend_right .= $db_subproblem;
            $db_seed .= $db_subproblem;
        }
        // extension left sequence
        else {
            $this->query_seq_extend_left = $query_subproblem . $this->query_seq_extend_left;
            $query_seed = $query_subproblem . $query_seed;
            $this->db_seq_extend_left = $db_subproblem . $this->db_seq_extend_left;
            $db_seed = $db_subproblem . $db_seed;
        }

        // use needleman algorithm to get optimal alignment
        $needlemanAlgorithm->get_optimal_alignment($query_seed, $db_seed);
        $needleman_optimal_alignment = $needlemanAlgorithm->get_needleman_opt_alignments();

        // get the max score for left and right extensions
        $this->optimal_match_score = $needleman_optimal_alignment['match_score'];
        $isOpt_max_score = $this->optimal_match_score >= max($this->score_array);
        if (empty($this->score_array) || $isOpt_max_score) {
            // extension right score
            if ($flag == FALSE) {
                $this->query_score_extend_right = $this->query_seq_extend_right;
                $this->db_score_extend_right = $this->db_seq_extend_right;
            }
            // extension left score
            else {
                $this->query_score_extend_left = $this->query_seq_extend_left;
                $this->db_score_extend_left = $this->db_seq_extend_left;
            }

            // implode returns a string from the elements of an array
            // implode ( array $pieces ) : string
            $query_needleman_alignment = implode($needleman_optimal_alignment['query_sequence']);
            $db_needleman_alignment = implode($needleman_optimal_alignment['$db_sequence']);
            $match_needleman_alignment = $needleman_optimal_alignment;
        }

        array_push($this->score_array, $this->optimal_match_score);

        $optimal_alignments = array(
            'query_seed' => $query_seed,
            'db_seed' => $db_seed,
            'query_needleman_alignment' => $query_needleman_alignment,
            'db_needleman_alignment' => $db_needleman_alignment,
            'match_needleman_alignment' => $match_needleman_alignment);

        return $optimal_alignments;
    }

    // set optimal alignments to data field
    private function set_optimal($optimal_alignments) {
        $this->query_seed = $optimal_alignments['query_seed'];
        $this->db_seed = $optimal_alignments['db_seed'];
        $this->query_needleman_alignment = $optimal_alignments['query_needleman_alignment'];
        $this->db_needleman_alignment = $optimal_alignments['db_needleman_alignment'];
        $this->match_needleman_alignment = $optimal_alignments['match_needleman_alignment'];
        return;
    }

}

// end class
?>

