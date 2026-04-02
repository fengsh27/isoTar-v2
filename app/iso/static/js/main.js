
	function checkVariableExistance(variable)
	{
		if (typeof variable === 'undefined' || variable === null) 
			return false;
		else return true;
	}
	/*
	$.validator.setDefaults({
    	errorElement: "span",
    	errorClass: "help-block",
    	highlight: function (element, errorClass, validClass) {
        	$(element).closest('.form-group').addClass('has-error');
    	},
    	unhighlight: function (element, errorClass, validClass) {
       		$(element).closest('.form-group').removeClass('has-error');
    	},
    	errorPlacement: function (error, element) {
    		console.log(element)
    		console.log(element.parent().length)
        	if (element.parent().length || element.prop('type') === 'checkbox' || element.prop('type') === 'radio') {
            	error.insertAfter(element.parent().parent());
        	} else {
            	error.insertAfter(element);
        	}
    	}
	});
	*/