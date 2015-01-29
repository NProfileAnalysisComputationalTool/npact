angular.module('npact')
  .config(function($provide) {
    $provide.decorator('$log', function($delegate, $sniffer) {
      var original = _.clone($delegate);
      var messages = [];
      _.forEach(['debug', 'error', 'info', 'log', 'warn'], function(lvl) {
        $delegate[lvl] = function() {
          messages.push(_.toArray(arguments));
          original[lvl].apply($delegate, arguments);
        };
      });
      $delegate.messages = messages;
      return $delegate;
    });

    $provide.decorator('$exceptionHandler', function($delegate, $injector) {
      return function(exception, cause) {
        var MessageBus = $injector.get('MessageBus');
        MessageBus.danger(exception + '(from "' + cause + '")');
        return $delegate(exception, cause);
    };

    });
  })

  .directive('npactMsgPane', function(MessageBus) {
    return {
      retrict: 'E',
      template: '<div id="msgpane" class="ui-corner-all"><npact-msg ng-repeat="msg in messages" message="msg"></npact-msg></div>',
      controller: function($scope, $log) {
        this.remove = _.bind(MessageBus.remove, MessageBus);
        $scope.$watch(
          function() { return MessageBus.messages; },
          function(val) {
            $scope.messages = val;
          });
      }
    };
  })
  .directive('npactMsg', function(STATIC_BASE_URL) {
    return {
      templateUrl: STATIC_BASE_URL + 'js/messages/message.html',
      require: '^npactMsgPane',
      scope: {message: '='},
      link: function($scope, $element, $attrs, ctrl) {
        $scope.dismiss = function() {
          $element.css({display: 'none'});
          ctrl.remove($scope.message);
        };
      }
    };
  })

  .service('MessageBus', function( $q, $timeout, $log) {
    'use strict';
    this.messages = [];

    this.remove = function(message) {
      this.messages = _.without(this.messages, message);
    };

    this.log = function(level, msg, hideWhen) {
      $log.log(level, msg);
      var obj = {level: level, text: msg, hideWhen: hideWhen};
      this.messages.push(obj);
      var remove = _.bind(this.remove, this, obj);
      if(hideWhen){
        if(isNaN(hideWhen)) {
          $q.when(hideWhen).then(remove);
        }
        else {
          $timeout(remove, hideWhen);
        }
      }
    };

    _.forEach(['info', 'danger', 'warning', 'success'], function(lvl) {
      this[lvl] = _.partial(this.log, lvl);
    }, this);
  })

  .service('EmailBuilder', function(MessageBus, $log, $location) {
    this.send = function() {
      var to = 'nathan@acceleration.net';
      var subject = 'NPACT run report';
      var url = $location.absUrl();
      var body = 'On page: \n' + url + '\n\n' + JSON.stringify($log.messages);
      var total = 'mailto:' + encodeURIComponent(to) +
            '?subject=' + encodeURIComponent(subject) +
            '&body=' + encodeURIComponent(body);
      window.location = total;
    };
  })
  ;
