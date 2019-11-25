<?php include './view/header.php'; ?>
<br>
<h1>Protein Subject: <?php echo trim($protein, ".txt"); ?></h1>
<h2>Sequences producing significant alignments:</h2>

<br>
<?php foreach ($alignments as $alignment) : ?>

    <label>
        Score:  
        <?php echo $alignment['optimal']['match_score']; ?>
    </label>
    <br>

    <table class="align">
        <tr>
            <td>Query </td>
            <td><?php echo $alignment['query_start_index'] + 1; ?></td>
            <td>
                <?php foreach ($alignment['optimal']['query_sequence'] as $v) : ?>
                <td>
                    <?php echo $v; ?>
                </td>
            <?php endforeach; ?>
            </td>
            <td> <?php echo $alignment['query_end_index'] + 1; ?></td>
        </tr>
        <tr>
            <td></td>
            <td></td>
            <td>
                <?php foreach ($alignment['optimal']['alignment_match'] as $v) : ?>
                <td>
                    <?php echo $v; ?>
                </td>
            <?php endforeach; ?>
            </td>
            <td> </td>
        </tr>
        <tr>
            <td>Subject </td>
            <td><?php echo $alignment['db_start_index'] + 1; ?></td>
            <td>
                <?php foreach ($alignment['optimal']['$db_sequence'] as $v) : ?>
                <td>
                    <?php echo $v; ?>
                </td>
            <?php endforeach; ?>
            </td>
            <td> <?php echo $alignment['db_end_index'] + 1; ?></td>
        </tr>
    </table>
    <br>
    <br>
<?php endforeach; ?>  

<?php include './view/footer.php'; ?>