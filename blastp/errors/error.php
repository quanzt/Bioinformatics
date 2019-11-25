<!DOCTYPE html>
<html>
    <body>
        <main>
            <!-- replacing separate database_error.php:
            report on error, with banner Database or class name if it's
                an Exception -->     
            <?php
            if ($alignment instanceof PDOException) { // including subclasses
                $label = 'Database';
                $error_message = $alignment->getMessage();
            } else if ($alignment instanceof Exception) {
                $label = get_class($alignment);
                $error_message = $alignment->getMessage();
            } else if (is_string($alignment)) {
                $label = 'User';
                $error_message = $alignment;
            } else {
                $label = 'Unclassified';
                $error_message = 'Error not Exception or string: bug in includer or error.php';
            }
            ?>
            <h1> <?php echo $label; ?> Error</h1>
            <h3> <?php echo $error_message; ?> </h3>
            <?php
            if ($alignment instanceof Exception) {
                echo '<p> at line ' . $alignment->getLine() . ' in file ' . $alignment->getFile() . '</p>';
                echo '<p> full backtrace:<br>';
                echo str_replace("\n", '<br>', $alignment->getTraceAsString());
                echo '</p>';
            }
            ?> 
            <p> <?php debug_print_backtrace(); ?></p>
        </main>
    </body>
</html>
