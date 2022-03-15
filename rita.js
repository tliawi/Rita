/*
rita.js

Copyright 2022 John R C Fairfield, see MIT License
    
Download, and put the two files index.html and rita.js in a folder (directory) and
open index.html in the Chrome browser.
*/

// rita ("ripple tank") creates a 2d wave propagation symplectic cellular automaton
// of the same size as the canvas, 
// and returns an object of methods to manipulate and display the cells

function rita(imgCanvas) {

var imgContext = imgCanvas.getContext("2d");
imgContext.globalAlpha = 1; //no transparency 
var height = imgCanvas.height;
var width = imgCanvas.width;
var maxi = height-1;
var maxj = width-1;
    
var imgData = imgContext.getImageData(0, 0, width, height); //is a copy
var rgba = imgData.data; //can modify this copy, and subsequently display with an imgContext.putImageData
var zoom = 1; //scale of drawImage
// Two dimensional structures (rgba and cellauto) are indexed in row major order (i.e. row, column) 
// with the origin in the upper left corner, 
// positive going down (in rows) and to the right (in cols).

var stiffness = 0.15;   
    
function setStiffness(b){
    stiffness = b;
}

function getParms(){
    return { stiffness:stiffness, clockTic: clockTic }
}

function allocate2d(v){
    var y = new Array(height);
    for (let i=0;i<height;i++) y[i] = (new Array(width)).fill(v);
    return y;
}
    
function clear(){
    clockTic = 0;
    for (let i=0;i<height;i++) for (let j=0;j<width;j++) u[i][j] = 0;
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
// from green (positive) to black (zero) to red (negative).
// Call with no third parm to have it adapt color range to datapoint having maximum magnitude.
// redBlack == 0 for red, 1 for Black
function renderToRGBA(numArr, redBlack, mag){ 
    
    function colorRGBA(num, ix) {
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
    
    if (mag == undefined){
        // find greatest magnitude, to normalize with
        mag = 0;
        for (let i=1;i<maxi;i++) for (let j=1; j<maxj; j++) if ( ((i+j)&1) == redBlack )
            mag = Math.max(mag,magnitude(numArr[i]));
    }
    else mag = Math.abs(mag); //guard against negative parameter
    if (mag==0) mag = 1; //avoid zerodivide in toxic all-zero case
    
    //rgba is linear canvas rgba data array of integer 0-255 values
    for (let i=1;i<maxi;i++) {
        for (let j=1; j<maxj; j++){
          let ix = 4*(i*width + j); //offset to R component in linear RGBA array
          if ( ((i+j)&1) == redBlack ){
            colorRGBA(numArr[i][j], ix );
          } else { //fill in with average of 4 neighbors of the same color
            colorRGBA( (numArr[i-1][j] + numArr[i+1][j] + numArr[i][j-1] + numArr[i][j+1])/4, ix);
          }
        }
    }
    
    //imgContext.putImageData(imgData,0,0,0,0,width*zoom,height*zoom); doesn't work.
    imgContext.putImageData(imgData, 0, 0); //any way to avoid doing both of these?
    imgContext.drawImage(imgCanvas,0,0,width*zoom,height*zoom);
    
}

var u = allocate2d(0.0);
    
function advanceT(){
    clockTic+= 2;
}

function getU()  {return u  }


    
// John Hartwell's W(t+1, L) = [ W(t, L) + a * W(t+1, L-1)]  / (1 + a), 
// where a = V * delta_t / delta_x
// Here Uij is current value of an absorbing cell, 
// and UPneighbor is an already computed up value of an interior cell
function absorbingFunc(Uij,UPneighbor){
    let a = Math.sqrt(beta/mass);
    return (Uij + a*UPneighbor)/(1+a);
}

// John Hartwell's second attempt, moving from t-0.5 to t,
// where a is still V * delta_t / delta_x. 
// W(t+1,L) = { 4*W(t,L) - W(t-1,L)  + a * [4*W(t+1,L-1) - W(t+1,L-2)] } / [3 * (1 + a)]
function absorbingFunc2(Uij, UMij, UPnbr, UPnbrnbr){
    let a = Math.sqrt(beta/mass);
    return (4*Uij - UMij + a*(4*UPnbr - UPnbrnbr))/(3*(1+a));
}


function computeRedBlack(){
    
    //red cells
    for (let i=1;i<maxi;i++) for (let j=1; j<maxj; j++) if ( !((i+j)&1) )  //note !
        u[i][j] += stiffness*(  (u[i+1][j]-u[i-1][j])  +  (u[i][j+1] - u[i][j-1])  ); //note +=
    
    //the red pass must be complete before the black pass, at time clockTic+1, is done
    
    //black cells
    for (let i=1;i<maxi;i++) for (let j=1; j<maxj; j++) if (  ((i+j)&1) )  //note no !
        u[i][j] -= stiffness*(  (u[i+1][j]-u[i-1][j])  +  (u[i][j+1] - u[i][j-1])  ); //note -=
    
}
    
function step(n=1){
    for (let i=0;i<n;i++){
        advanceT();
        computeRedBlack();
    }
}
    
function setZoom(z=1){ zoom = z; }

// d is distance from center, sigma > 0 is variance, width, of bump
function gaussian(distanceSquared,sigma){
    const c = Math.sqrt(2*Math.PI);
    return Math.exp(-(distanceSquared/(sigma*sigma))/2)/(sigma*c);
}

//Per Ostrov and Rucker is bad to seed with discontinuities
//Does not seed edges. Seeds red or black cells in u, which after advanceT will be u and um
function seedGaussian(row,col,redBlack,sigma,amplitude=1){
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++) if ( ((i+j)&1) == redBlack) { //does not work if parens around lhs omitted!
            let g = amplitude*gaussian((i-row)*(i-row)+(j-col)*(j-col), sigma);
            u[i][j]  += g; 
            //continuous in space and time, 
            //in that is at top dead center, about to turn down, 
            //i.e. vertical velocity is zero
        }
    }
}

/*
//does not cover borders
function laplace(ux,i,j){
    //NW, SE, NE, SW
    return ux[i-1][j-1] + ux[i+1][j+1] + ux[i-1][j+1] + ux[i+1][j-1] - 4*ux[i][j];
}
    
function secondDerivativeLaplace(uIn, uOut){
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            uOut[i][j] = laplace(uIn,i,j);
        }
    }
}

//seeds with second derivative of gaussian.
//Does not seed edges. Seeds up and u, which after advanceT will be u and um
//cannot be used during animation, since it clobbers up rather than adds to it
function seedMexicanHat(row,col,redBlack,sigma,amplitude=1){
    seedGaussian(row,col,redBlack,sigma,amplitude); //adds gaussian hill to up and u
    secondDerivativeLaplace(u,up,1); //clobbers up to laplacian of u
    copyU(up,u); //copy laplacian back into u
}
*/
    
function integral(ux){
    let r = 0;
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            r += ux[i][j];
        }
    }
    return r;
}

var clockTic = 0; //counts step advanceT's
function setClockTic(st){clockTic=st}

/*
var pumpCol=0;
function setPumpCol(pc){pumpCol=pc; setClockTic(-1);}
var pumpPeriod=9;
function setPumpPeriod(p){pumpPeriod=p; setClockTic(-1);}
    
function pumpColumn(someU){ //period is number of tics in full cycle
    if (pumpCol<=0 || pumpCol>=width-1) return; //can't pump borders
    let amp = Math.sin(2*Math.PI*clockTic/pumpPeriod);
    for (let i=1;i<height-1;i++){ someU[i][pumpCol] = amp;}
}
*/
    
return { 
        clear:clear,
        clearRGBA: clearRGBA,
        imgContext:imgContext,
        imgData:imgData,
        renderToRGBA: renderToRGBA,
        setZoom:setZoom,
        setStiffness:setStiffness,
        getParms:getParms,
        step: step,
        seedGaussian: seedGaussian,
       // seedMexicanHat: seedMexicanHat,
        integral: integral,
       
        setClockTic: setClockTic,
        //setPumpCol: setPumpCol,
        //setPumpPeriod: setPumpPeriod,
        getU:getU, 
    
    };
}
