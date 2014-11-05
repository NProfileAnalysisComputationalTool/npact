angular.module('npact', ['infinite-scroll'])
  .value('K', Kinetic)
  .value('Err', {
    TrackAlreadyDefined: new Error('Track with this name already defined, each track name must be unique'),
    TrackNotFound: new Error('track not found'),
    ProfileNotFound: new Error('profile not found')
  })
  .value('Evt', {
    REDRAW:'redraw',
    NOOP:'noop',
    REBUILD:'rebuild',
    PAN:'pan',
    ZOOM:'zoom'
  })

  .run(function($rootScope, STATIC_BASE_URL) {
    $rootScope.STATIC_BASE_URL = STATIC_BASE_URL;
  })
;
