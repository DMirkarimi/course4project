//Setting variable up for parallax effect
const woop = document.getElementById("test");

//Decreasing scrollspeed for parallax effect

window.addEventListener("scroll", function()
    {
    let offset = window.pageYOffset;
    woop.style.backgroundPositionY = offset / 10 + "px";
    })