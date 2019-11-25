<?php include './view/header.php'; ?>

<html>
    <main>
        <section>
            <form autocomplete="off" action="index.php" method="post">
                <input type="hidden" name="action" value="getBlastp_result">
                <h1> BLAST: blastp suite </h1>
                <div class="tab">
                    <h2>Standard Protein BLAST</h2>
                    <br>
                    <label>Query:</label><br>
                    <p class="tab">
                        <input type='text' name="query_seq" required="required">
                    </p>
                    <br>
                    <br>
                    <label>Protein Subject:</label><br>
                    <div class="tab">
                        <select name="protein">
                            <option value="4AUE_A.txt">4AUE_A</option>
                            <option value="2PTN_A.txt">2PTN_A</option>                     
                            <option value="4B31_A.txt">4B31_A</option>
                        </select>
                    </div>
                </div>
                <br><br>
                <h1> Algorithm parameters </h1>
                <div class="tab">
                    <h2>General Parameters</h2>
                    <div class="tab">
                        <label>Score threshold:</label>
                        <input type='number' name="score_threshold" value="10" style="width: 4em" required="required" >
                        (Stop extension if alignment score fall below score threshold.)
                        <br>
                        <br>
                        <label>Word size:</label>
                        <input type='number' name="word_size" value="6" style="width: 4em" required="required" >
                        (The length of the seed that initiates an alignment.)
                    </div>
                    <br><br>

                    <h2>Scoring Parameters</h2>
                    <div class="tab">
                        <label>Matrix:</label><br>
                        <div class="tab">
                            <select name="matrix">
                                <option value="blosum62">BLOSUM62</option>
                                <option value="blosum50">BLOSUM50</option>
                                <option value="pam250">PAM250</option>                     
                            </select>
                        </div>
                        (Assigns a score for aligning pairs of residues, and determines overall alignment score.)
                        <br><br>
                        <label>Gap Costs:</label>
                        <input type='number' name="gap_costs"  value="-1" style="width: 4em" required="required" >
                        (Cost to create and extend a gap in an alignment.)
                    </div>
                </div>
                <br><br>
                <input type="submit" 
                       style="height: 50px; width: 150px; left: 250; top: 250; font-face: 'Comic Sans MS'; font-size: larger; color: teal; background-color: #FFFFC0; border: 3pt ridge lightgrey"
                       value="BLAST" />
                <br><br>
            </form>
        </section>
    </main>
</html>

<?php include './view/footer.php'; ?>