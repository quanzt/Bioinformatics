<?php

// Reference: some modification from https://github.com/aebruno/needleman-wunsch/blob/master/needleman-wunsch-class.php
class NeedlemanWunsch {

    private static $arrow_up = '&#8593;';
    private static $arrow_left = '&#8592;';
    private static $arrow_nw = '&#8598;';
    private $alignment_matrix = array();
    private $alignments = array();
    private $matrix = array();
    private $gap_cost = -1;

    public function __construct($gap_cost, $matrix) {
        $this->matrix = $matrix;
        $this->gap_cost = $gap_cost;
    }

    // generate the global aligment 
    public function get_optimal_alignment($query_alignment, $db_alignment) {
        $this->setup_matrix_grid($query_alignment, $db_alignment);
        $scoring_matrix = $this->scoring_system($this->matrix);
        for ($i = 1; $i < count($this->alignment_matrix); $i++) {
            for ($j = 1; $j < count($this->alignment_matrix[$i]); $j++) {
                //table score e.g. BLOSUM62
                $match_score = $scoring_matrix[self::convert_char_int($query_alignment[$i - 1])][self::convert_char_int($db_alignment[$j - 1])];
                $match = $this->alignment_matrix[$i - 1][$j - 1]['val'] + $match_score;
                $hgap = $this->alignment_matrix[$i - 1][$j]['val'] + $this->gap_cost;
                $vgap = $this->alignment_matrix[$i][$j - 1]['val'] + $this->gap_cost;
                $max = max($match, $hgap, $vgap);
                $pointer = self::$arrow_nw;
                if ($max === $hgap) {
                    $pointer = self::$arrow_up;
                } else if ($max === $vgap) {
                    $pointer = self::$arrow_left;
                }
                $this->alignment_matrix[$i][$j]['pointer'] = $pointer;
                $this->alignment_matrix[$i][$j]['val'] = $max;
            }
        }
        $i = count($this->alignment_matrix) - 1;
        $j = count($this->alignment_matrix[0]) - 1;
        $this->alignments['query_sequence'] = array();
        $this->alignments['$db_sequence'] = array();
        $this->alignments['alignment_match'] = array();
        $this->alignments['match_score'] = $this->alignment_matrix[$i][$j]['val'];
        while ($i !== 0 and $j !== 0) {
            $query_sequence = $query_alignment[$i - 1];
            $db_sequence = $db_alignment[$j - 1];
            $this->alignment_matrix[$i][$j]['trace'] = true;
            $pointer = $this->alignment_matrix[$i][$j]['pointer'];
            if ($pointer === self::$arrow_nw) {
                $i--;
                $j--;
                $this->alignments['query_sequence'][] = $query_sequence;
                $this->alignments['$db_sequence'][] = $db_sequence;
                $this->alignments['alignment_match'][] = ($query_sequence === $db_sequence) ? $query_sequence : ' ';
            } else if ($pointer === self::$arrow_up) {
                $i--;
                $this->alignments['query_sequence'][] = $query_sequence;
                $this->alignments['$db_sequence'][] = '-';
                $this->alignments['alignment_match'][] = ' ';
            } else if ($pointer === self::$arrow_left) {
                $j--;
                $this->alignments['query_sequence'][] = '-';
                $this->alignments['$db_sequence'][] = $db_sequence;
                $this->alignments['alignment_match'][] = ' ';
            } else {
                die("Invalid pointer: $i,$j");
            }
        }
        foreach (array('query_sequence', '$db_sequence', 'alignment_match') as $k) {
            $this->alignments[$k] = array_reverse($this->alignments[$k]);
        }
        return $this->alignment_matrix;
    }

    // Constructing the grid
    private function setup_matrix_grid($query_alignment, $db_alignment) {
        $this->alignment_matrix = array();
        $this->alignments = array();
        for ($i = 0; $i < strlen($query_alignment) + 1; $i++) {
            for ($j = 0; $j < strlen($db_alignment) + 1; $j++) {
                $this->alignment_matrix[$i][$j] = array(
                    'pointer' => null,
                    'trace' => null,
                    'val' => 0
                );
            }
        }
        for ($i = 0; $i < strlen($query_alignment); $i++) {
            $this->alignment_matrix[$i + 1][0]['val'] = ($i + 1) * $this->gap_cost;
        }
        for ($j = 0; $j < strlen($db_alignment); $j++) {
            $this->alignment_matrix[0][$j + 1]['val'] = ($j + 1) * $this->gap_cost;
        }
    }

    // Choosing a scoring system and fill out the matrix
    public function scoring_system($substituion_matrices) {
        $matrix = array();
        $maxtrix_scoring = file('./models/' . $substituion_matrices . '.txt');
        foreach ($maxtrix_scoring as $score_fill_out) {
            $matrix[] = explode(',', $score_fill_out);
        }
        return $matrix;
    }

    // getter method
    public function get_needleman_opt_alignments() {
        return $this->alignments;
    }

    // Protein sequences are chains of letters drawn from a 22-letter 
    // alphabet of amino acids (20 standard plus 2 unusual ones
    public static function convert_char_int($char) {
        if ($char == 'A') {
            return 0;
        } else if ($char == 'R') {
            return 1;
        } else if ($char == 'N') {
            return 2;
        } else if ($char == 'D') {
            return 3;
        } else if ($char == 'C') {
            return 4;
        } else if ($char == 'Q') {
            return 5;
        } else if ($char == 'F') {
            return 6;
        } else if ($char == 'G') {
            return 7;
        } else if ($char == 'H') {
            return 8;
        } else if ($char == 'E') {
            return 9;
        } else if ($char == 'L') {
            return 10;
        } else if ($char == 'K') {
            return 11;
        } else if ($char == 'M') {
            return 12;
        } else if ($char == 'T') {
            return 13;
        } else if ($char == 'P') {
            return 14;
        } else if ($char == 'S') {
            return 15;
        } else if ($char == 'V') {
            return 16;
        } else if ($char == 'W') {
            return 17;
        } else if ($char == 'Y') {
            return 18;
        } else if ($char == 'I') {
            return 19;
        } else {
            echo "Not An Amino Acid Alphabet";
            return;
        }
    }

}

?>