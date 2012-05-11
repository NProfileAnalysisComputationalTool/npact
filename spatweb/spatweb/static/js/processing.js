//set up spindowns	

$.fn.spinDown = function() {
	return this.click(
            function() {
		var $this = $(this);
		$this.next().slideToggle();
		$this.find('.ui-icon').toggleClass('ui-icon-triangle-1-s');
		return false;
		
	});
	
};

// The .bind method from Prototype.js 
if (!Function.prototype.bind) { // check if native implementation available
  Function.prototype.bind = function(){ 
    var fn = this, args = Array.prototype.slice.call(arguments),
        object = args.shift(); 
    return function(){ 
      return fn.apply(object, 
        args.concat(Array.prototype.slice.call(arguments))); 
    }; 
  };
}


var progress = {
    start: null,
    
    init: function(url) {
        progress.url = url;
        progress.startRequest();
        progress.timerDom = jQuery('#timer > span');
        progress.steps =  jQuery('#steps');
        $('#progress  h3').spinDown();
    },

    updateStatusDisplay: function(data) {
        //{"steps": ['a','b'], "tdiff": 10.32341, "step_desc": "a step description."}
        console.log('Update Progress Display', data);
        if(!data) return;
        var steps = data.steps;
        for(var i=0; i < steps.length; i++) {
            progress.steps.prepend('<li>' + data.steps[i] + '</li>');
        }
        var files = data.files;
        if (files && files.length) {
            var newnode = $('<ul id="files"/>');
            for(i=0; i < files.length; i++){
                var href = files[i];
                var name = /([^/]*)$/.exec(href)[1];
                $(newnode).append('<li><a href="' + href + '">' + name + '</a></li>');
            }
            $('#files').replaceWith(newnode);
        }
    },

    startRequest: function(){
        console.log('Starting Progress Request Chain');
        if(!progress.start) {
            progress.start = new Date();
            progress.interval = window.setTimeout(progress.updateTimer, 1000);
        }
        jQuery.ajax({ url: progress.url,
                      dataType:'json',
                      success: progress.onsuccess,
                      error: progress.onerror
                    });
    },

    onerror: function(data){
        console.log('Error', data);
        if(progress.interval)
            clearTimeout(progress.interval);

        setTimeout(function() { progress.updateStatusDisplay(jQuery.parseJSON(data.responseText)); }, 0);
        $('#computationrunning').slideToggle()
            .after(function () {
                       $('#computationerred').slideToggle();
                   });
    },

    onsuccess: function(data){
        console.log('Progress', data);
        //async call
        setTimeout(progress.updateStatusDisplay.bind(progress, data));

        if(data.next=='results')
            progress.processResults(data);
        else 
            progress.startRequest();

    },
    
    processResults: function(data) {
        console.log('finished: ', data.download_url);
        if(progress.interval)
            clearTimeout(progress.interval);
        $('#title').html('N-PACT Is Finished');
        $('#downloadlink').attr('href', data.download_url);
        $('#timer').toggleClass('ui-state-highlight');
        $('#progressreport').trigger('click');
        $('#results').fadeIn(100);
        $('#results h3').trigger('click');
        $('#computationrunning').fadeOut(100);
        //setTimeout(function() {window.location = data.download_url; }, 800);
    },
    
    updateTimer: function (){
        if(progress.start) {
            $('#timer').fadeIn();
            var elapsed = new Date() - progress.start;
            progress.timerDom.html(Math.floor(elapsed/1000) + ' sec.');
            progress.interval = setTimeout(progress.updateTimer, 1010 - (elapsed % 1000));
        }
    }
};
