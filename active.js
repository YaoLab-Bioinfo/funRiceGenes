
$(document).ready(function() {

    /* clear file button control */
    var fileControl = $("#file1");

    $("#clearFile1").on("click", function () {
        fileControl.replaceWith( fileControl = fileControl.clone( true ) );
        $("#file1_progress").hide();
    });

    $("#uploadFormat").on("change", function () {
        fileControl.replaceWith( fileControl = fileControl.clone( true ) );
        $("#file1_progress").hide();
    });

    /* file input progress bar control */
    $( "#file1" ).change(function() {
      document.getElementById("file1_progress").setAttribute('style', "height:20px; margin-top:5px;");
    });

});


