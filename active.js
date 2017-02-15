
$(document).ready(function() {

    /* clear file button control */
    var fileControl = $("#genfamin");

    $("#clear6").on("click", function () {
        fileControl.replaceWith( fileControl = fileControl.clone( true ) );
        $("#genfamin_progress").hide();
    });

    /* file input progress bar control */
    $( "#genfamin" ).change(function() {
      document.getElementById("genfamin_progress").setAttribute('style', "height:20px; margin-top:5px;");
    });

});


