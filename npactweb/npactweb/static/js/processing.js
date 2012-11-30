//set up spindowns

$.fn.spinDown = function() {
    return this.click(
            function() {
        var $this = $(this);
        $this.next().slideToggle();
        //$this.find('.ui-icon').toggleClass('ui-icon-triangle-1-s');
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
    POLLTIME: 2500,

    start: null,

    init: function(url) {
        progress.url = url;
        window.setTimeout(progress.startRequest, 50);
        progress.currentstepdesc = $('#currentstepdesc');
    },

    updateStatusDisplay: function(data) {
        //{"steps": ['a','b'], "tdiff": 10.32341, "step_desc": "a step description."}
        console.log('Update Progress Display', data);
        if(!data) return;

        var steps = data.steps;
        if(steps) {
            progress.currentstepdesc.html(steps.join("<br/>"));
        }
        // for(var i=0; i < steps.length; i++) {
        //     progress.steps.prepend('<li>' + data.steps[i] + '</li>');
        // }

        var files = data.files;
        if (files && files.length) {
            var newnode = $('<ul id="files"/>');
            for(i=0; i < files.length; i++){
                var href = files[i];
                var name = /([^/]*)$/.exec(href)[1];
                $(newnode).append('<li><a href="' + href + '" target="_blank">' + name + '</a></li>');
            }
            $('#files').replaceWith(newnode);
            console.log("Finished building file list.");
            $('#results').fadeIn(200);
        }
    },

    startRequest: function(){
        console.log('Starting Progress Request Chain');
        if(!progress.start) {
            progress.start = new Date();

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
                $('.onerror').slideToggle();
            });
    },

    onsuccess: function(data){
        console.log('Progress', data);
        //async call
        if(data.next=='results')
            setTimeout(progress.processResults.bind(progress, data));
        else
            setTimeout(progress.startRequest.bind(progress, data), progress.POLLTIME);

        setTimeout(progress.updateStatusDisplay.bind(progress, data));
    },

    processResults: function(data) {
        console.log('finished: ', data.download_url);
        if(progress.interval)
            clearTimeout(progress.interval);
        $('#downloadlink').attr('href', data.download_url);
        $('#title').html('N-PACT Is Finished');
        $('#progressreport').fadeOut(100);
        $('.onsuccess').fadeIn(600);

    },
};
