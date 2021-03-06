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
    
var alpha = 0; // vertical spring to 0, alpha*u .  0.002?
var beta = 0.2;    // horizontal springs to 4-neighbors, laPlacian beta*(N+S+E+W - 4u)
var mass = 1;
    
function setAlpha(a){
    alpha = a;
}
function setBeta(b){
    beta = b;
}
    
function setMass(m){
    mass = m;
}

function getParms(){
    return { alpha:alpha, beta:beta, mass:mass , clockTic: clockTic }
}

function allocate2d(v){
    var y = new Array(height);
    for (let i=0;i<height;i++) y[i] = (new Array(width)).fill(v);
    return y;
}
    
function clear(){
    potentialSides = 0;
    potentialMain  = 0;
    clockTic = -1;
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
// from green (positive) to black (zero) to red (negative).
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

var cellEnergy = allocate2d(0.0);
    
//does not cover border cells
//assumes all three of up, u, and um (at clockTic+1, clockTic, and clockTic-1) are valid
function getEnergies(){
    
    let potentialSides = 0;
    let potentialMain  = 0;
    let kinetic = 0;
    
    var emitted = 0;
    var absorbed = 0;
    
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            let Uij = u[i][j];
            let del1 = up[i][j] - Uij;
            let del2 = Uij - um[i][j]; 
            
            let N = (u[i-1][j] - Uij);
            let S = (u[i+1][j] - Uij);
            let E = (u[i][j+1] - Uij);
            let W = (u[i][j-1] - Uij);
            let NSEW = N+S+E+W; //laplacian u[i-1][j] + u[i+1][j] + u[i][j+1] + u[i][j-1] - 4*Uij;
            
            //divide by two because you're calculating 
            //work done against the spring to take it to current position, 
            //and spring's force increases linearly with displacement,
            //so taking triangular half of square area force*force,
            //since choosing units such that 
            //spatial (i,j) units equal u (orthogonal spring displacement) units
            
            let epMain = alpha*Uij*Uij/2; //potential energy in main spring displacement
            // each spring beta*y*y/2 is counted twice, by two neighbors, so divide by 2 again
            let epSides = beta*(N*N+S*S+E*E+W*W)/4; //in side springs displacements
            
            //This is magic. Not del1^2, not del2^2, but del1*del2. This can go negative, 
            // at top and bottom of swing, when del changes directions
            // and dels are smallest in magnitude (so is a very small negative value).
            //This magic makes total energy work out.
            //See second sheet of bit.ly/testAccessOpen simple harmonic oscillator spreadsheet
            let eKinetic = mass*del1*del2/2;
            
            potentialMain += epMain;
            potentialSides += epSides;
            kinetic += eKinetic;
            
            cellEnergy[i][j] = epMain + epSides + eKinetic;
            
        }
    }
    
    return { 
        potentialSides:potentialSides,
        potentialMain:potentialMain,
        kinetic:kinetic,
        emitted:emitted,
        absorbed:absorbed,
        cellEnergy:cellEnergy,
    }
}

var fudge=1;
function setFudge(f){fudge = f;}

// John Hartwell's W(t+1, L) = [ W(t, L) + a * W(t+1, L-1)]  / (1 + a), 
// where a = V * delta_t / delta_x
// Here Uij is current value of an absorbing cell, 
// and UPneighbor is an already computed up value of an interior cell
function absorbingFunc(Uij,UPneighbor){
    let a = Math.sqrt(beta/mass)*fudge;
    return (Uij + a*UPneighbor)/(1+a);
}

// John Hartwell's second attempt, moving from t-0.5 to t,
// where a is still V * delta_t / delta_x. 
// W(t+1,L) = { 4*W(t,L) - W(t-1,L)  + a * [4*W(t+1,L-1) - W(t+1,L-2)] } / [3 * (1 + a)]
function absorbingFunc2(Uij, UMij, UPnbr, UPnbrnbr){
    let a = Math.sqrt(beta/mass)*fudge;
    return (4*Uij - UMij + a*(4*UPnbr - UPnbrnbr))/(3*(1+a));
}

var funcSwitch=false;
//true or false, absorbingFunc or absorbingFunc2
function toggleFuncSwitch(tf){funcSwitch = !funcSwitch; return funcSwitch;} 
    
// calculates up = u + du/dt - alpha*u + beta*laplacian
function standardFunc(Uij,Umij,N,S,E,W){ 

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
            + alpha*(-Uij)/mass // main spring accelleration
            //+ beta*(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1] -4*u[i][j]) // Laplacian
            // N,S,E,W are *differences* betwen u and its 4 neighbors
            + beta*(N+S+E+W)/mass; //side springs can cancel each other, do so on any flat slope
}

function computeUp(){
    
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
    
    /*
    //four edges
    //left and right
    for (let i=0;i<maxi;i++) { 
        up[i][0] = noBorder(i,0);
        up[i][maxj] = noBorder(i,maxj);
    }
    //top and bottom
    for (let j=0;j<maxj;j++){
        up[0][j] = noBorder(0,j);
        up[maxi][j] = noBorder(maxi,j);
    }
    */
    
    
    //core interior cells, not on border
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            let UIJ = u[i][j];
            up[i][j] = func[i][j](UIJ, um[i][j],
                       u[i-1][j]-UIJ, u[i+1][j]-UIJ, u[i][j+1]-UIJ, u[i][j-1]-UIJ);
        }
    }
    
    
    //four edges
    //for absorbingFuncs, must be done AFTER interior cells have been UPped
    if (funcSwitch){
        
        //left and right
        for (let i=1;i<maxi;i++) { 
            up[i][0] = absorbingFunc(u[i][0], up[i][1]);
            up[i][maxj] = absorbingFunc(u[i][maxj], up[i][maxj-1]);
        }
        //top and bottom
        for (let j=1;j<maxj;j++){
            up[0][j] = absorbingFunc(u[0][j], up[1][j]);
            up[maxi][j] = absorbingFunc(u[maxi][j], up[maxi-1][j]);
        }
        //kludge corners
        up[0][0] = absorbingFunc(u[0][0], up[1][1]);
        up[maxi][0] = absorbingFunc(u[maxi][0], up[maxi-1][1]);
        up[0][maxj] = absorbingFunc(u[0][maxj], up[1][maxj-1]);
        up[maxi][maxj] = absorbingFunc(u[maxi][maxj], up[maxi-1][maxj-1]);
        
    } else {
        
        //left and right
        for (let i=1;i<maxi;i++) { 
            up[i][0] = absorbingFunc2(u[i][0], um[i][0], up[i][1], up[i][2]);
            up[i][maxj] = absorbingFunc2(u[i][maxj], um[i][maxj], up[i][maxj-1], up[i][maxj-2]);
        }
        //top and bottom
        for (let j=1;j<maxj;j++){
            up[0][j] = absorbingFunc2(u[0][j], um[0][j], up[1][j], up[2][j]);
            up[maxi][j] = absorbingFunc2(u[maxi][j], um[maxi][j], up[maxi-1][j], up[maxi-2][j]);
        }
        //kludge corners
        up[0][0] = absorbingFunc2(u[0][0], um[0][0], up[1][1], up[2][2]);
        up[maxi][0] = absorbingFunc2(u[maxi][0], um[maxi][0], up[maxi-1][1], up[maxi-2][1]);
        up[0][maxj] = absorbingFunc2(u[0][maxj], um[0][maxj], up[1][maxj-1], up[1][maxj-2]);
        up[maxi][maxj] = absorbingFunc2(u[maxi][maxj], um[maxi][maxj], up[maxi-1][maxj-1], up[maxi-2][maxj-2]);
        
    }

    pumpColumn(up); //does nothing if pumpCol == 0, else overwrites up[][pumpCol].
    
}  
    
function step(n=1){
    for (let i=0;i<n;i++){
        advanceT(); //um <= u; u <= up; clockTic++
        computeUp();
        //so now up, u and um are valid, needed for getEnergies() of u
        //renderToRGBA(u); // when we display u, we've already computed up as well
    }
}
    
function setZoom(z=1){ zoom = z; }

// d is distance from center, sigma > 0 is variance, width, of bump
function gaussian(distanceSquared,sigma){
    const c = Math.sqrt(2*Math.PI);
    return Math.exp(-(distanceSquared/(sigma*sigma))/2)/(sigma*c);
}

//Per Ostrov and Rucker is bad to seed with discontinuities
//Does not seed edges. Seeds up and u, which after advanceT will be u and um
function seedGaussian(row,col,sigma,amplitude=1){
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            let g = amplitude*gaussian((i-row)*(i-row)+(j-col)*(j-col), sigma);
            u[i][j]  += g; 
            up[i][j] += g;
            //continuous in space and time, 
            //in that is at top dead center, about to turn down, 
            //i.e. vertical velocity is zero
        }
    }
}

//does not cover borders
//one use of amplifier is -1
function copyU(uIn,uOut, amplifier=1){
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            uOut[i][j] = amplifier*uIn[i][j] ;
        }
    }
}

//does not cover borders
function laplace(ux,i,j){
    //N, S, E, W
    return ux[i-1][j] + ux[i+1][j] + ux[i][j+1] + ux[i][j-1] - 4*ux[i][j];
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
function seedMexicanHat(row,col,sigma,amplitude=1){
    seedGaussian(row,col,sigma,amplitude); //adds gaussian hill to up and u
    secondDerivativeLaplace(u,up,1); //clobbers up to laplacian of u
    copyU(up,u); //copy laplacian back into u
}

function integral(ux){
    let r = 0;
    for (let i=1;i<maxi;i++) { 
        for (let j=1;j<maxj;j++){
            r += ux[i][j];
        }
    }
    return r;
}

var clockTic = -1; //counts step advanceT's
function setClockTic(st){clockTic=st}
    
var pumpCol=0;
function setPumpCol(pc){pumpCol=pc; setClockTic(-1);}
var pumpPeriod=9;
function setPumpPeriod(p){pumpPeriod=p; setClockTic(-1);}
    
function pumpColumn(someU){ //period is number of tics in full cycle
    if (pumpCol<=0 || pumpCol>=width-1) return; //can't pump borders
    let amp = Math.sin(2*Math.PI*clockTic/pumpPeriod);
    for (let i=1;i<height-1;i++){ someU[i][pumpCol] = amp;}
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
        setMass: setMass,
        getParms:getParms,
        getEnergies: getEnergies,
        step: step,
        seedGaussian: seedGaussian,
        seedMexicanHat: seedMexicanHat,
        integral: integral,
        standardFunc: standardFunc,
        setClockTic: setClockTic,
        setPumpCol: setPumpCol,
        setPumpPeriod: setPumpPeriod,
        getUp: getUp, getU:getU, getUm: getUm,
    
        setFudge:setFudge,
        absorbingFunc:absorbingFunc,
        toggleFuncSwitch:toggleFuncSwitch,
        absorbingFunc2:absorbingFunc2,
    };
}
