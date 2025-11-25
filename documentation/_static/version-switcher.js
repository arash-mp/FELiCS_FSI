fetch('/versions.json')
  .then(r => r.json())
  .then(versions => {
    const select = document.createElement('select');
    versions.forEach(v => {
      const option = document.createElement('option');
      option.value = v.path;
      option.text = v.name;
      if (window.location.pathname.startsWith(v.path) && v.path !== '/') option.selected = true;
      select.appendChild(option);
    });
    select.onchange = () => { window.location.href = select.value; };
    document.querySelector('header')?.appendChild(select);  // adjust for your theme
  });
