(function () {
  'use strict';

  var VERSION_SEGMENT_PATTERN = /^v[0-9]+\.[0-9]+(\.[0-9]+)?$/;

  function findVersionRoot(pathname) {
    var segments = pathname.split('/').filter(Boolean);
    var versionIndex = segments.findIndex(function (segment) {
      return segment === 'development' || VERSION_SEGMENT_PATTERN.test(segment);
    });
    if (versionIndex === -1) {
      return null;
    }
    var baseSegments = segments.slice(0, versionIndex);
    return {
      currentVersion: segments[versionIndex],
      // No leading "/" when there is no prefix: a bare pagesBase of "/" would
      // turn "pagesBase + '/switcher.json'" into "//switcher.json", which
      // browsers parse as a protocol-relative URL (host "switcher.json")
      // instead of an absolute path.
      pagesBase: baseSegments.length ? '/' + baseSegments.join('/') : '',
    };
  }

  function buildSwitcherBar(versions, currentVersion) {
    var bar = document.createElement('div');
    bar.className = 'felics-version-switcher';

    var label = document.createElement('span');
    label.className = 'felics-version-switcher__label';
    label.textContent = 'Version:';
    bar.appendChild(label);

    var select = document.createElement('select');
    select.className = 'felics-version-switcher__select';
    versions.forEach(function (entry) {
      var option = document.createElement('option');
      option.value = entry.url;
      option.textContent = entry.name;
      option.selected = entry.version === currentVersion;
      select.appendChild(option);
    });
    select.addEventListener('change', function () {
      window.location.href = select.value;
    });
    bar.appendChild(select);

    return bar;
  }

  function init() {
    var versionRoot = findVersionRoot(window.location.pathname);
    if (versionRoot === null) {
      return;
    }
    fetch(versionRoot.pagesBase + '/switcher.json')
      .then(function (response) {
        if (!response.ok) {
          throw new Error('switcher.json request failed: ' + response.status);
        }
        return response.json();
      })
      .then(function (versions) {
        var bar = buildSwitcherBar(versions, versionRoot.currentVersion);
        document.body.appendChild(bar);
      })
      .catch(function (error) {
        console.warn('FELiCS version switcher could not be loaded:', error);
      });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
