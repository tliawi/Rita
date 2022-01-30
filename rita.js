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
var maxi = height-1;
var maxj = width-1;
    
var imgData = imgContext.getImageData(0, 0, width, height); //is a copy
var rgba = imgData.data; //can modify this copy, and subsequently display with an imgContext.putImageData
var zoom = 1; //scale of drawImage
// Two dimensional structures (rgba and cellauto) are indexed in row major order (i.e. row, column) 
// with the origin in the upper left corner, 
// positive going down (in rows) and to the right (in cols).

// C+ = gain*((N0+S0+E0+W0) + present*C0 + past*C- )
    
var alpha = 0.002; // vertical spring to 0, alpha*u 
var beta = 0.2;    // horizontal springs to 4-neighbors, laPlacian beta*(N+S+E+W - 4u)

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

var um = allocate2d(0.0); //u minus, prior time step
var u  = allocate2d(0.0); //current
var up = allocate2d(0.0); //u plus, next time step, calculated as fcn of u and um
    
var func = allocate2d(standardFunc); //determines the function to be applied at that pixel
    
function advanceT(){
    var z = um;
    um = u;
    u  = up;
    up = z; //to be overwritten subsequently
    clockTic++;
}
    
function getUp() {return up }
function getU()  {return u  }
function getUm() {return um }

//global energies
var potentialSides = 0;
var potentialMain  = 0;
var emitted = 0;
var absorbed = 0;
    
function getEnergies(){
    return { 
        potentialSides:potentialSides,
        potentialMain:potentialMain,
        kinetic:kineticEnergy(),
        emitted:emitted,
        absorbed:absorbed,
    }
}

function kineticEnergy(){
    let e = 0;
    for (let i=0;i<height;i++) for (let j=0;j<width;j++) {
        let del = u[i][j]-um[i][j];
        e+= del * del/2; //not *beta! Spring tension's got nothing to do with kinetic energy
    }
    return e;
}
    
// calculates up = u + du/dt - alpha*u + beta*laplacian
function standardFunc(Uij,Umij,N,S,E,W){ 

    //update energies
    //divide by two because you're calculating 
    //work done against the spring to take it to current position, 
    //and spring's force increases linearly with displacement,
    //so taking triangular half of square area force*force,
    //since choosing units such that 
    //spatial (i,j) units equal u (orthogonal spring displacement) units
    potentialMain  += alpha*Uij*Uij/2; //in main spring displacement
    potentialSides += beta*(N*N + S*S + E*E + W*W)/2; //in side springs displacements

    //first term is u
    //remaining terms will make up next iteration's du/dt
    //second term is current du/dt
    //alpha and beta reflect spring stiffness/mass
    //alpha term is accel (change in du/dt) from main, central vertical spring at i,j
    //beta term is accel (change in du/dt) from springs connecting i,j to NSEW neighbors,
    //must be zero on flat slope 
    return  Uij //current value
            + Uij - Umij  //current du/dt
            // plus changes in du/dt for next iteration:
            + alpha*(-Uij) // main spring displacement force
            //+ beta*(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1] -4*u[i][j]) // Laplacian
            // N,S,E,W are *differences* betwen u and its 4 neighbors
            + beta*(N+S+E+W); //side springs can cancel each other, do so on any flat slope
}

function computeUp(){
    
    potentialSides = 0;
    potentialMain  = 0;

    //At boundries met by annulling counterwave?
    //Suppose 'missing' 4-neighbors at boundaries are virtualized as
    //the negative of their opposing side (thinking NSEW, if N is missing, use -S, etc).
    //But that turns out to be the same as fixing the entire boundary at 0, since 4 corners will be 0
    //(N-S, E-W, so 0), and 4 edges will be independent of what happens in interior--for example, a N-S edge
    //will be insensitive to anything in interior (since E-W is zero), so will only propagate along the edge,
    //which if it wasn't initialized to nonzero values will stay 0.
    //returns u[i][j] save for border cells, when it returns inner neighbor
    
    function noBorder(i,j){
        i = i==0?1:i; 
        i = i==maxi?maxi-1:i;
        
        j = j==0?1:j;
        j = j==maxj?maxj-1:j;
        
        return u[i][j];
    }
    
    //four edges
    //left and right
    for (let i=0;i<maxi;i++) { 
        up[i][0] = noBorder(i,0);
        up[i][maxj] = noBorder(i,maxj);
    }
    //top and bottom
    for (let j=1;j<maxj;j++){
        up[0][j] = noBorder(0,j);
        up[maxi][j] = noBorder(maxi,j);
    }
    
    //core
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            let UIJ = u[i][j];
            up[i][j] = func[i][j](UIJ, um[i][j],
                       u[i-1][j]-UIJ, u[i+1][j]-UIJ, u[i][j+1]-UIJ, u[i][j-1]-UIJ);
        }
    }
    
    pumpColumn(); //does nothing if pumpCol == 0, else overwrites up[][pumpCol].
    
}  
    
function step(n=1){
    for (let i=0;i<n;i++){
        computeUp();
        advanceT();
        renderToRGBA(u);
    }
}
    
function setZoom(z=1){ zoom = z; }

// d is distance from center, sigma > 0 is variance, width, of bump
function gaussian(distanceSquared,sigma){
    const c = Math.sqrt(2*Math.PI);
    return Math.exp(-(distanceSquared/(sigma*sigma))/2)/(sigma*c);
}

//Per Ostrov and Rucker is bad to seed with discontinuities
//Does not seed edges
function seed(row,col,sigma,amplitude=1){
    var maxi = height - 1;
    var maxj = width - 1;
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            um[i][j] = u[i][j] = amplitude*gaussian((i-row)*(i-row)+(j-col)*(j-col), sigma); //continuous in space
            //and time, in that is at top dead center, about to turn down, i.e. vertical velocity is zero
        }
    }
}

var clockTic=0; //counts step advanceT's
function setClockTic(st){clockTic=st}
var pumpCol=0;
function setPumpCol(pc){pumpCol=pc; setClockTic(0);}
var pumpPeriod=9;
function setPumpPeriod(p){pumpPeriod=p; setClockTic(0);}
    
function pumpColumn(){ //period is number of tics in full cycle
    if (pumpCol<=0 || pumpCol>=width-1) return; //can't pump borders
    let amp = Math.sin(2*Math.PI*clockTic/pumpPeriod);
    for (let i=0;i<height-1;i++){ up[i][pumpCol] = amp;}
}

    
return { 
        clear:clear,
        clearRGBA: clearRGBA,
        imgContext:imgContext,
        imgData:imgData,
        renderToRGBA: renderToRGBA,
        setZoom:setZoom,
        setAlpha: setAlpha,
        setBeta: setBeta,
        getEnergies: getEnergies,
        kineticEnergy: kineticEnergy,
        step: step,
        seed: seed,
        setClockTic: setClockTic,
        setPumpCol: setPumpCol,
        setPumpPeriod: setPumpPeriod,
        getUp: getUp, getU:getU, getUm: getUm,
    };
}
