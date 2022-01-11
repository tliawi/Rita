/*
rita.js

Copyright 2022 John R C Fairfield, see MIT License
    
Download, and put the two files index.html and rita.js in a folder (directory) and
open index.html in the Chrome browser.
*/

// rita ("ripple tank") creates a 2d wave propagation cellular automaton of the same size as the canvas, 
// and returns an object of methods to manipulate and display the cells

function rita(imgCanvas) {

var imgContext = imgCanvas.getContext("2d");
imgContext.globalAlpha = 1; //no transparency  
var imgData = imgContext.getImageData(0, 0, imgCanvas.width, imgCanvas.height); //is a copy
var rgba = imgData.data; //can modify this copy, and subsequently display with an imgContext.putImageData
    
// Two dimensional structures (rgba and cellauto) are indexed in row major order (i.e. row, column) 
// with the origin in the upper left corner, 
// positive going down (in rows) and to the right (in cols).
    

function allocate2d(){
    var y = new Array(imgContext.height);
    for (let i=0;i<imgContext.height;i++) y[i] = (new Array(imgContext.width)).fill(0);
    return y;
}

function clearRGBA(){
    for (let i=0;i<imgCanvas.height;i++) for (let j=0;imgCanvas.width;j++)  {
            let ix = 4*(i*imgCanvas.width + j);
            rgba[ix] = 0;
            rgba[ix+1] = 0;
            rgba[ix+2] = 0;
            rgba[ix+3] = 0xff;
    }
    imgContext.putImageData(imgData, 0, 0);
}

// clobbers rgba, a global r,g,b,alpha linear canvas data array visualizing the numeric data 
// in a given complete two dimensional numeric Array
// Numeric values are rendered in a spectrum 
// from green (positive) to grey (zero) to red (negative).
function renderToRGBA(numArr, threshold=9999999999){
    const GREYVAL = 64;
    var height = numArr.length, width = numArr[0].length;
    
    // find greatest magnitude, to normalize with
    var mag = 0;
    for (let i=0;i<height;i++) mag = Math.max(mag,magnitude(numArr[i]));
    if (mag==0) mag = 1; //avoid zerodivide in toxic all-zero case
    //console.log("magnitude "+mag);
    
    //normalize very differently if threshold is so small that unit differences will be clearly distinguishable.
    if (threshold <= 10) mag = threshold;
    
    //rgba is linear canvas rgba data array of integer 0-255 values
    
    for (let i=0;i<height;i++){
        for (let j=0;j<width;j++){
            let ix = 4*(i*width + j);
            let num = numArr[i][j];
            if (num<0) {
                num = num < -threshold? -threshold: num;
                rgba[ix]   = Math.floor(((255-GREYVAL)*(-num))/mag) + GREYVAL; //red negative
                rgba[ix+1] =  GREYVAL;   //green
            } else {
                num = num > threshold? threshold: num;
                rgba[ix]   =  GREYVAL;   //red
                rgba[ix+1] = Math.floor(((255-GREYVAL)*num)/mag) + GREYVAL; //green positive
            }
            rgba[ix+2] = GREYVAL; //blue
            rgba[ix+3] = 255; //alpha
        }
    }
    
    imgContext.putImageData(imgData, 0, 0);
}

um = allocate2d(); //u minus, prior time step
u  = allocate2d(); //current
up = allocate2d(); //u plus, next time step, calculated as fcn of u and um
    
meta = allocate2d(); //determines the function to be applied at that pixel
    
//current velocity = (u - um)/T, where T is time it takes to do iteration, here taken as 1

function advanceT(){
    var z = um;
    um = u;
    u  = up;
    up = z; //to be overwritten subsequently
}


function computeUp(){
    maxi = u.length - 1;
    maxj = u[0].length - 1;
    
    //four corners
    up[0][0]       = (u[1][0]+u[0][1]) - um[0][0];
    up[maxi][maxj] = (u[maxi][maxj-1]+u[maxi-1][maxj]) - um[maxi][maxj];
    up[0][maxj]    = (u[0][maxj-1]+u[1][maxj]) - um[0,maxj];
    up[maxi][0]    = (u[maxi-1][0]+u[maxi][1]) - um[maxi][0];
    
    //four edges
    for (int i=1;i<maxi;i++) { 
        up[i][0] = (u[i-1][0]+u[i+1][0]+u[i,1]) - 2*um[i][0];
        up[i][maxj] = (u[i-1][maxj]+u[i+1][maxj]+u[i][maxj-1]) - 2*um[i][maxj];
    }
    for (int j=1;j<maxj;j++){
        up[0][j] = (u[0][j-1]+u[0][j+1]+u[1,j]) - 2*um[0][j];
        up[maxi][j] = (u[maxi][j-1]+u[maxi][j+1]+u[maxi-1][j]) - 2*um[maxi][j];
    }
    
    //core
    for (int i=1;i<maxi;i++) for (int j=1;j<maxj;j++){
        up[i][j] = (u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1]) - 3*um[i][j];
    }
    
}  
    
function step(n=1){
    for (let i=0;i<n;i++){
        computeUp();
        advanceT();
        renderToRGBA(u);
    }
}
    
return { 
        clearImage: clearImage,
        clearRGBA: clearRGBA,
        renderToRGBA: renderToRGBA,
        step: step,
        um: um,
        u: u,
        up:up,
    };
}
