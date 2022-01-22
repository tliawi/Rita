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
var height = imgCanvas.height;
var width = imgCanvas.width;
    
var imgData = imgContext.getImageData(0, 0, width, height); //is a copy
var rgba = imgData.data; //can modify this copy, and subsequently display with an imgContext.putImageData
var zoom = 1; //scale of drawImage
// Two dimensional structures (rgba and cellauto) are indexed in row major order (i.e. row, column) 
// with the origin in the upper left corner, 
// positive going down (in rows) and to the right (in cols).

// C+ = gain*((N0+S0+E0+W0) + present*C0 + past*C- )
    
var alpha = 0.01;
var beta = 0.02;

function setAlpha(a){
    alpha = a;
}
function setBeta(b){
    beta = b;
}

function allocate2d(v){
    var y = new Array(height);
    for (let i=0;i<height;i++) y[i] = (new Array(width)).fill(v);
    return y;
}
    
function clear(){
    potentialSides = 0;
    potentialMain  = 0;
    for (let i=0;i<height;i++) for (let j=0;j<width;j++) u[i][j] = up[i][j] = um[i][j] = 0;
    clearRGBA();
}

function clearRGBA(){
    for (let i=0;i<height;i++) for (let j=0;j<width;j++)  {
            let ix = 4*(i*width + j);
            rgba[ix] = 0;
            rgba[ix+1] = 0;
            rgba[ix+2] = 0;
            rgba[ix+3] = 0xff;
    }
    imgContext.putImageData(imgData, 0, 0);
}

    
function magnitude(numArr){ //numArr a one dimensional numeric array
    return Math.max(Math.max(...numArr), -Math.min(...numArr));
}

// clobbers rgba, a global r,g,b,alpha linear canvas data array visualizing the numeric data 
// in a given height x width numeric array (like um,u,up)
// Numeric values are rendered in a spectrum 
// from green (positive) to yellow (zero) to red (negative).
// Call with no second parm to have it adapt color range to datapoint having maximum magnitude.
function renderToRGBA(numArr, mag){
    
    if (mag == undefined){
        // find greatest magnitude, to normalize with
        mag = 0;
        for (let i=0;i<height;i++) mag = Math.max(mag,magnitude(numArr[i]));
    }
    else mag = Math.abs(mag); //guard against negative parameter
    if (mag==0) mag = 1; //avoid zerodivide in toxic all-zero case
    
    //rgba is linear canvas rgba data array of integer 0-255 values
    for (let i=0;i<height;i++){
        for (let j=0;j<width;j++){
            let ix = 4*(i*width + j);
            let num = numArr[i][j];
            if (num<0) {
                num = num < -mag? mag: -num; //num now positive <= mag
                let s = Math.floor(255*num/mag); //runs from 0 to 255 as num from 0 to mag
                rgba[ix]   = s; //red
                rgba[ix+1] = Math.floor(s/2); //green
                rgba[ix+2] = s; //blue
            } else {
                num = num > mag? mag: num; //num positive <= mag
                let s = Math.floor(255*num/mag);//runs from 0 to 255 as num from 0 to mag
                rgba[ix]   = Math.floor(s/2); //red
                rgba[ix+1] = s; //green
                rgba[ix+2] = s; //blue
            }
            rgba[ix+3] = 255; //alpha
        }
    }
    
    //imgContext.putImageData(imgData,0,0,0,0,width*zoom,height*zoom); doesn't work.
    imgContext.putImageData(imgData, 0, 0); //any way to avoid doing both of these?
    imgContext.drawImage(imgCanvas,0,0,width*zoom,height*zoom);
    
}

um = allocate2d(0.0); //u minus, prior time step
u  = allocate2d(0.0); //current
up = allocate2d(0.0); //u plus, next time step, calculated as fcn of u and um
    
//meta = allocate2d(0); //determines the function to be applied at that pixel
    
//current velocity = (u - um)/T, where T is time it takes to do iteration, here taken as 1

function advanceT(){
    var z = um;
    um = u;
    u  = up;
    up = z; //to be overwritten subsequently
}
    
function getUp() {return up }
function getU()  {return u  }
function getUm() {return um }

var potentialSides = 0;
var potentialMain  = 0;
    
function getPotentialSides(){
    return potentialSides;
}

function getPotentialMain(){
    return potentialMain;
}
    
function kineticEnergy(){
    var e = 0;
    for (let i=0;i<height;i++) for (let j=0;j<width;j++) {
        let del = u[i][j]-um[i][j];
        e+= del * del/2; //not *beta!
    }
    return e;
}

function computeUp(){
    var maxi = height - 1;
    var maxj = width - 1;
    
    potentialSides = 0;
    potentialMain  = 0;
    
    //average of my neighbors current values, minus my prior value
    
    //at boundries met by annulling counterwave.
    //Suppose 'missing' 4-neighbors at boundaries are virtualized as
    //the negative of their opposing side (thinking NSEW, if N is missing, use -S, etc).
    //But that turns out to be the same as fixing the entire boundary at 0, since 4 corners will be 0
    //(N-S, E-W, so 0), and 4 edges will be independent of what happens in interior--for example, a N-S edge
    //will be insensitive to anything in interior (since E-W is zero), so will only propagate along the edge,
    //which if it wasn't initialized to nonzero values will stay 0.
    
    /*
    //four corners
    up[0][0]       = (u[1][0]+u[0][1])/2.0                  - um[0][0];
    up[maxi][maxj] = (u[maxi][maxj-1]+u[maxi-1][maxj])/2.0  - um[maxi][maxj];
    up[0][maxj]    = (u[0][maxj-1]+u[1][maxj])/2.0          - um[0][maxj];
    up[maxi][0]    = (u[maxi-1][0]+u[maxi][1])/2.0          - um[maxi][0];
    
    //four edges
    for (let i=1;i<maxi;i++) { 
        up[i][0] = (u[i-1][0]+u[i+1][0]+u[i][1])/3.0                - um[i][0];
        up[i][maxj] = (u[i-1][maxj]+u[i+1][maxj]+u[i][maxj-1])/3.0  - um[i][maxj];
    }
    for (let j=1;j<maxj;j++){
        up[0][j] = (u[0][j-1]+u[0][j+1]+u[1][j])/3.0                - um[0][j];
        up[maxi][j] = (u[maxi][j-1]+u[maxi][j+1]+u[maxi-1][j])/3.0  - um[maxi][j];
    }
    */
    //four corners
    up[0][0]       = 0;
    up[maxi][maxj] = 0;
    up[0][maxj]    = 0;
    up[maxi][0]    = 0;
    
    //four edges
    for (let i=1;i<maxi;i++) { 
        up[i][0] = (u[i-1][0]+u[i+1][0]+u[i][1])/3.0                - um[i][0];
        up[i][maxj] = (u[i-1][maxj]+u[i+1][maxj]+u[i][maxj-1])/3.0  - um[i][maxj];
    }
    for (let j=1;j<maxj;j++){
        up[0][j] = (u[0][j-1]+u[0][j+1]+u[1][j])/3.0                - um[0][j];
        up[maxi][j] = (u[maxi][j-1]+u[maxi][j+1]+u[maxi-1][j])/3.0  - um[maxi][j];
    }
    
    //core
    // up = u + du/dt + alpha*laplacian - beta*u
    var N,S,W,E,UIJ;
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            //first term is u
            //remaining terms will make up next iteration's du/dt
            //second term is current du/dt
            //alpha term is accel (change in du/dt) from Laplacian springs connecting i,j to NSEW neighbors, 
            //must be zero on flat slope 
            //beta term is accel (change in du/dt) from main, central vertical spring at i,j
            //alpha and beta reflect spring stiffness/mass
            UIJ = u[i][j];
            N = u[i-1][j]-UIJ; 
            S = u[i+1][j]-UIJ; 
            W = u[i][j-1]-UIJ; 
            E = u[i][j+1]-UIJ;
            potentialSides += alpha*(N*N + S*S + E*E + W*W)/2; //in side springs displacements
            potentialMain  += beta*UIJ*UIJ/2; //in main spring displacement
            up[i][j] =  UIJ 
                        + UIJ - um[i][j]  //current du/dt
                        //+ alpha*(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1] -4*u[i][j]) //alpha * Laplacian
                        + alpha*(N+S+E+W) //side springs can cancel each other, do so if flat slope
                        - beta*UIJ; //beta * main spring displacement
        }
    }
    
}  
    
function step(n=1){
    for (let i=0;i<n;i++){
        computeUp();
        advanceT();
        renderToRGBA(u);
    }
}
    
function setZoom(z=1){ zoom = z; }
    
return { 
        clear:clear,
        clearRGBA: clearRGBA,
        imgContext:imgContext,
        imgData:imgData,
        renderToRGBA: renderToRGBA,
        setZoom:setZoom,
        setAlpha: setAlpha,
        setBeta: setBeta,
        getPotentialSides: getPotentialSides,
        getPotentialMain: getPotentialMain,
        kineticEnergy: kineticEnergy,
        step: step,
        getUp: getUp, getU:getU, getUm: getUm,
    };
}
