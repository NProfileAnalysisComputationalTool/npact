angular.module('npact', ['ngMessages', 'sticky', 'ngSanitize', 'dialogService', 'ngCookies'])
  .config(function($locationProvider) {
    'use strict';
    $locationProvider.html5Mode(true);
  })
  .run(function($rootScope, STATIC_BASE_URL) {
    'use strict';
    $rootScope.STATIC_BASE_URL = STATIC_BASE_URL;
  })

  .directive('contenteditable', ['$sce', function($sce) {
    'use strict';
    return {
      restrict: 'A', // only activate on element attribute
      require: '?ngModel', // get a hold of NgModelController
      link: function(scope, element, attrs, ngModel) {
        if (!ngModel) return; // do nothing if no ng-model

        // Specify how UI should be updated
        ngModel.$render = function() {
          element.html($sce.getTrustedHtml(ngModel.$viewValue || ''));
        };

        // Listen for change events to enable binding
        element.on('blur keyup change', function() {
          scope.$evalAsync(read);
        });
        read(); // initialize

        // Write data to the model
        function read() {
          var html = element.html();
          // When we clear the content editable the browser leaves a <br> behind
          // If strip-br attribute is provided then we strip this out
          if (attrs.stripBr) {
            html = html.replace(/<br>$/, '');
          }
          ngModel.$setViewValue(html);
        }

        scope.$on('$destroy', function() {
          element.off('blur keyup change');
        });
      }
    };
  }])

  .directive('checkboxList', function() {
    'use strict';
    return {
      restrict: 'E',
      scope: { options: '=' },
      require: 'ngModel',
      template: '<label ng-repeat="key in options">{{key}}<input type="checkbox" ng-model="selected[key]"></label>',
      link: function(scope, element, attrs, ngModel) {
        scope.selected = _.object(scope.options);
        ngModel.$render = function() {
          //start with everything false
          scope.selected = _.object(scope.options);
          //set selected list to true
          _.forEach(ngModel.$viewValue, function(el) {
            scope.selected[el] = true;
          });
        }.bind(this);

        element.on('change', function() {
          scope.$evalAsync(function() {
            var newVal = [];
            _.forEach(scope.selected, function(v, k) {
              if(v) {
                newVal.push(k);
              }
            });
            ngModel.$setViewValue(newVal);
          });
        });
        scope.$on('$destroy', function() { element.off('change'); });
      }
    };
  })

  .directive('qtipHelpIcon', function($log) {
    'use strict';
    return {
      restrict: 'C',
      link: function($scope, $element, $attrs) {
        $log.log("In link function: ", $element.next('.qtip-help-text'));
        $element.qtip({
           content: {
             text: $element.parent().next('.qtip-help-text')
            }
        });
      }
    };
  })
;
