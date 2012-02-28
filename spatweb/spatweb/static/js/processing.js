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

var progress = {
    start: null,
    
    init: function(url) {

        progress.url = url;
        progress.startRequest();
        progress.timerDom = jQuery('#timer > span');
        progress.steps =  jQuery('#steps');
        $('#progress  h3').spinDown();
    },

    updateDisplay: function(data) {
        //{"steps": ['a','b'], "tdiff": 10.32341, "step_desc": "a step description."}
        console.log('Update Progress Display', data);
        if(!data) return;
        for(var i=0; i< data.steps.length; i++) {
            progress.steps.append('<li>' + data.steps[i] + '</li>');
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
        if(progress.interval)
            window.clearInterval(progress.interval);

        $('#computationrunning').slideToggle()
            .after(function () {
                       $('#computationerred').slideToggle();
                   });
    },

    onsuccess: function(data){
        console.log('Progress', data);
        if(data.next=='results'){
            progress.processResults(data);
        }
        else {
            //data.pt has progress data
            progress.startRequest();
            progress.updateDisplay(data.pt);
        }
    },
    
    processResults: function(data) {
        console.log('finished: ', data.download_url);
        if(progress.interval)
            clearTimeout(progress.interval);
        $('#title').html('N-PACT Is Finished');
        $('#downloadlink').attr('href', data.download_url);
        $('#reconfigurelink').attr('href', data.reconfigure_url);
        $('#progressreport').trigger('click');
        $('#results').fadeIn(100);
        $('#results h3').trigger('click');
        $('#computationrunning').fadeOut(100);
        setTimeout(function() {window.location = data.download_url; }, 800);
    },
    
    updateTimer: function (){
        if(progress.start) {
            $('#timer').fadeIn();
            var elapsed = new Date() - progress.start;
            progress.timerDom.html(Math.floor(elapsed/1000) + ' sec.');
            progress.interval = setTimeout(progress.updateTimer, elapsed % 1000);
        }
    }
};
