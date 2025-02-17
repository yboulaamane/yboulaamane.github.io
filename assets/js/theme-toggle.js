document.addEventListener("DOMContentLoaded", function () {
  const toggleSwitch = document.querySelector(".theme-switch input");
  const body = document.body;

  // Load user preference
  if (localStorage.getItem("theme") === "dark") {
    body.classList.add("dark-mode");
    toggleSwitch.checked = true;
  }

  // Toggle dark mode
  toggleSwitch.addEventListener("change", function () {
    if (this.checked) {
      body.classList.add("dark-mode");
      localStorage.setItem("theme", "dark");
    } else {
      body.classList.remove("dark-mode");
      localStorage.setItem("theme", "light");
    }
  });
});
