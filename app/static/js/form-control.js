function populateParts(){

    var model_class = document.getElementById("model_class_id");
    var part_type_id = document.getElementById("part_type_id");
    var selected_class = model_class.options[model_class.selectedIndex].value;

    if (selected_class=='Intein-based logic gate model')
    {
        //metrics.options.length=0;
        part_type_id.options[1].disabled = true;
        part_type_id.options[2].disabled = false;
        part_type_id.options[3].disabled = true;
        part_type_id.options[4].disabled = true;
        part_type_id.options[5].disabled = true;
        part_type_id.options[6].disabled = true;
        part_type_id.options[7].disabled = true;
    }
    else if (selected_class=='ECF/AS-based logic gate model')
    {
        //metrics.options.length=0;
        part_type_id.options[1].disabled = true;
        part_type_id.options[2].disabled = true;
        part_type_id.options[3].disabled = true;
        part_type_id.options[4].disabled = false;
        part_type_id.options[5].disabled = true;
        part_type_id.options[6].disabled = true;
        part_type_id.options[7].disabled = true;
    }
    else if (selected_class=='Create custom model') {
        window.location="/custom.html";
    }
    else {
        part_type_id.options[1].disabled = false;
        part_type_id.options[2].disabled = false;
        part_type_id.options[3].disabled = false;
        part_type_id.options[4].disabled = false;
        part_type_id.options[5].disabled = false;
        part_type_id.options[6].disabled = false;
        part_type_id.options[7].disabled = false;
    }
}

function showGuide() {
    alert("Guidelines for formatting the input files:\n\n- Please save your file in a .csv format.\n- Please use the first-row for the variable/gate names.\n- Please use the first-column as the index (or timepoints in case of dynamic data).\n- Each n-consecutive rows will be treated as all induction states for an n-inputs logic gate.");
}

function validateEmail(email) {
    const re = /^(([^<>()[\]\\.,;:\s@"]+(\.[^<>()[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;
    return re.test(String(email).toLowerCase());
}

function validateFile(filename) {

    // Allowing file type
    var allowedExtensions = /(\.csv)$/i;
    //var allowedExtensions = /(\.csv|\.xlsx|\.xls)$/i;
    if (!allowedExtensions.exec(filename)) {
        return false;
    }
    return true;
}

function validateModelForm() {

    var msg = ""
    var err_count = 0

    var model = document.getElementById("model_class_id").selectedIndex;
    var part = document.getElementById("part_type_id").selectedIndex;
    var fluo = document.getElementById("fluorescence_data_id").value;
    var od = document.getElementById("od_data_id").value;
    var name = document.getElementById("model_name_id").value;
    
    if(model==0) {
        msg += "Please select a valid model!\n";
        err_count++;
    }
    if(part==0) {
        msg += "Please select a valid genetic part!\n";
        err_count++;
    }
    if(fluo=="") {
        msg += "Please provide fluorescence data!\n";
        err_count++;
    }
    if(od==""){
        msg += "Please provide cell-growth data!\n";
        err_count++;
    }
    if(fluo!="" && !validateFile(fluo)) {
        msg += "Make sure the fluorescence data format is .csv!\n";
        err_count++;
    }
    if(od!="" && !validateFile(od)) {
        msg += "Make sure the cell-growth data format is .csv!!\n";
        err_count++;
    }
    if(name=="") {
        msg += "Please provide the model name!!\n";
        err_count++;
    }

    if(err_count>0){
        alert(msg)
        event.preventDefault();
    }
    else {
        email = prompt("This task will take some time to run.\nPlease enter your email address, to receive the report once it is done:", "")
        if (email==null) {
            event.preventDefault();
        }
        else if (email=="" || !validateEmail(email)) {
            alert("Please enter a valid email address!");
            event.preventDefault();
        } else if (validateEmail(email)) {
            alert("Your task is submitted. Thank you for using Xcalibr!\nThe result is coming up soon...");
        }
    }
}

function validateCustomForm() {

    var msg = ""
    var err_count = 0

    var name = document.getElementById("model_name_id").value;
    var part = document.getElementById("part_type_id").selectedIndex;
    var crn = document.getElementById("crn_id").value;
    
    if(name=="") {
        msg += "Please provide the model name!!\n";
        err_count++;
    }
    if(part==0) {
        msg += "Please select a valid genetic part!\n";
        err_count++;
    }
    if(crn=="") {
        msg += "Please specify the reaction networks!\n";
        err_count++;
    }
    
    if(err_count>0){
        alert(msg)
        event.preventDefault();
    }
    else {
        alert("New model is added!");
    }
}

function validateDesignForm() {

    var msg = ""
    var err_count = 0

    var model = document.getElementById("model_id").selectedIndex;
    
    if(model==0 || model==(document.getElementById("model_id").length - 1)) {
        msg += "Please select a valid model!!\n";
        err_count++;
    }

    if(err_count>0){
        alert(msg)
        event.preventDefault();
    }
    else {
        email = prompt("This task will take some time to run.\nPlease enter your email address, to receive the report once it is done:", "")
        if (email==null) {
            event.preventDefault();
        }
        else if (email=="" || !validateEmail(email)) {
            alert("Please enter a valid email address!");
            event.preventDefault();
        } else if (validateEmail(email)) {
            alert("Your task is submitted. Thank you for using Xcalibr!\nThe result is coming up soon...");
        }
    }
}