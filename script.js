var out = null;
var ardl = null;
var isrun = false;
var likelihood={y: [], x:[], name: 'likelihood', type: 'scatter', yaxis: 'y2',  mode: 'lines'};
var residuals={y: [], x:[], name: 'residuals', type:'scatter', yaxis:'y1', mode:'lines'};
//var prettyprices={y: [], x:[], name: 'price', type:'scatter', yaxis:'y2', mode:'marker'};
var fifth={y: [], x:[], name: '5th quantile', type:'scatter', yaxis:'y1', mode:'lines'};
var ninetyfifth={y: [], x:[], name: '95th quantile', type:'scatter', yaxis:'y1', mode:'lines'};
var mean_={y:[],x:[], name:'mean', type:'scatter', yaxis:'y1',mode:'lines'};
var fitted = {y:[],x:[], name:'Modelled', type:'scatter', yaxis:'y2', mode:'lines'};
var q95 = {y:[],x:[], name:'95th Quantile of Residuals', type:'scatter', yaxis:'y2',mode:'lines'};
var q05 = {y:[],x:[], name:'5th Quantile of Residuals', type:'scatter', yaxis:'y2',mode:'lines'};
var q95UCI={y: [], x:[], name: '95th quantile upper 95% CI', type:'scatter', yaxis:'y1', mode:'lines'};
var q05LCI={y: [], x:[], name: '5th quantile lower 95% CI', type:'scatter', yaxis:'y1', mode:'lines'};
var q95UCIexp={y: [], x:[], name: '95th quantile upper 95% CI', type:'scatter', yaxis:'y2', mode:'lines'};
var q05LCIexp={y: [], x:[], name: '5th quantile lower 95% CI', type:'scatter', yaxis:'y2', mode:'lines'};
var prices={y: [], x:[], name: 'price', type:'scatter', yaxis:'y2', mode:'markers'};

    var rvf = {y:[],x:[], name:'Fitted v Residuals', type:'scatter', yaxis:'y',mode:'markers'};
    var ll, lr, pr, diff, meandistance, tomorrow;
var meanres ,stdev ,p05 ,p95,p05_L ,p95_U ;

var owa_baseUrl = 'https://btconometrics.com/owa/';
var owa_cmds = owa_cmds || [];
owa_cmds.push(['setSiteId', '90acb5e6721f2bf8f364cb0cb3aa7e6a']);
owa_cmds.push(['trackPageView']);
owa_cmds.push(['trackClicks']);

(function() {
	var _owa = document.createElement('script'); _owa.type = 'text/javascript'; _owa.async = true;
	owa_baseUrl = ('https:' == document.location.protocol ? window.owa_baseSecUrl || owa_baseUrl.replace(/http:/, 'https:') : owa_baseUrl );
	_owa.src = owa_baseUrl + 'modules/base/js/owa.tracker-combined-min.js';
	var _owa_s = document.getElementsByTagName('script')[0]; _owa_s.parentNode.insertBefore(_owa, _owa_s);
}());

function deepClone(X) {
  return JSON.parse(JSON.stringify(X));
}

function rescale(maxnew, minnew, maxold, minold, value) {
  return (((maxnew-minnew)/(maxold-minold))*(value-maxold)+maxnew);
}

/*
 * pre: 
 *  x is a vector 
 *  centre is a pivot point
 *
 * post:
 *  returns an array filled with the likelihood of each point in x, given the history of x
 */
function getLikelihood(x, centre) {
  var c = centre|0;
  var output =[];
  for(var i = 0; i< x.length; ++i) {
    var ov = x.slice(0,i).filter(counted=>counted>=x[i]).length;
    var uv = x.slice(0,i).filter(counted=>counted<=x[i]).length;
    if(x[i]<c) {
        output.push(getNPLikelihood(uv+1,i+1,0.05)); //+1 shift from zero for true df
      } else if(x[i]>=c) {
        output.push(getNPLikelihood(ov+1,i+1,0.05));
      }
    }
    return output;
}

function Wald(thetahats, theta0s, varthetahats) {
 var diff = math.subtract(thetahats, theta0s);
 var tp = math.transpose(diff )
 var varinv = math.multiply(tp, math.inv(varthetahats))
 var W = math.multiply(varinv, diff);
 return W;
  
}

function egranger(y, x) {
  
  var molslm = ARDL(y,x);
  var olslm = molslm.lm
  var residlm = olslm.residuals;
  var Dresidlm = diff(residlm,1);
  var Lresidlm = lag(residlm,1);//.slice(1,Dresid.length+1);
  var n = Dresidlm.length;
  var egcmlm = mOLS(Dresidlm,Lresidlm, false); //no intercept
  var taus=egcmlm.beta1_tstat;
  return taus;
}
  


//returns the [by] lag of x
function lag(_x,by) {
  var x = deepClone(_x);
  return x.slice(0,x.length-by);
}

function diff(_x, by) {
  var x = deepClone(_x);
  var d= [];
  for(i = by;i<x.length;i++) {
    d[i-by] = x[i]-x[i-by]
  }
  return d;
}

function minimum(value) {
      if (toString.call(value) !== "[object Array]")
      return false;
      return Math.min.apply(null, value);
   }

function maximum(value) {
      if (toString.call(value) !== "[object Array]")
      return false;
      return Math.max.apply(null, value);
   }
  //implements a simple percentile function
  //x is a vector, p is the probability
  function percentile(_x, p) {
    //because we dont want to be mutating x
    var x = deepClone(_x);
//take a shortcut if p >= 1 or p<= 0
  if(p>=1) {
      return (maximum(x));
    } if (p<=0) {
      return (minimum(x));
      }else {
    //first get the rank
    var n = x.length;
    var r = (1+ p*(n-1))-1; //-1 for zero addressing

    //now split the rank into integer and fractional components
    function f(n){ return Number("0."+String(n).split('.')[1] || 0); }
    function i(n){ return Number(String(n).split('.')[0] || 0); }

    var rfrac = f(r);
    var rint = i(r);

    //now sort the vector ascending
    x.sort(function(a, b){return a-b});

    //now interpolate for the correct position:
    var xp = (1-rfrac)*x[rint]+rfrac*x[rint+1];
    return xp;
  }
  }

function rmNaN(_X) {
  var x = deepClone(_X);
  xnan  = x.filter(function (value) {return !Number.isNaN(value);});
  return xnan
}

//stepping through the set of reals to find the solution
//probably could go to nlogn if i did binary search but meh
function binominv(n,p,a) {
  var setOfReals=[];
  //big enough for most purposes!
  for (j=0;j<1000000000000;j++) {
   var result = jStat.binomial.cdf(j,n,p);
      if(result>a) { 
        console.log(j);
        return(j);
      }
    }
}

// gives a 1-alpha upper CI on the pth quantile of X
/*
 * Lower quantile confidence interval
=percentile(X,(binom.inv (n,p,1-alpha/2)+1)/n

Upper quantile confidence interval
=percentile(X,( binom.inv (n,p,alpha/2)+1)/n


Where 
?	X is vector
?	n is length X
?	p is the quantile of interest
?	alpha is 1-coverage
?	binom.inv parameters = (trials, probability_success, alpha)
*/
function upperCI(alpha, p, _x) {
  var X = deepClone(_x);
  var n = X.length;
  return percentile(X,binominv(n,p,1-alpha/2)/n);
}

//taking the complement for the lower bound
function lowerCI(alpha, p, _x) {
  var X = deepClone(_x);
  var n = X.length;
  return percentile(X,binominv(n,p,alpha/2)/n);
}


  function mean(_x) {
   var X = deepClone(_x);
  return (X.reduce((previous, current) => current += previous)/X.length);
}

//returns the (non-parametric) likelihood that there are x nonconforming units in n
    function getNPLikelihood(x,n, alpha) {
      return jStat.beta.inv( alpha, x, n-x+1 );
    }
    function plot1(data, T="Residuals") {
      // console.log(data);
       var layout = {
         title: 'S2F BTC',
         yaxis: {title: T},
         xaxis: {
           type: 'date'
         }, 
         showlegend: true,
        legend: {orientation: 'h',
    bgcolor: 'transparent'}
       };
    Plotly.newPlot('chart', data, layout).then(postPlotHandler('loader_c1'));

     }
     function postPlotHandler(divid){
         document.getElementById(divid).style.display = "none";
       };

var quantdata="";

     function plot2(data2, title_, yaxistitle_, yaxis2title_, divid_) {
        var _title = title_ || 'BTC S2F'
        var _yaxistitle = yaxistitle_ || 'Residuals'
        var _yaxis2title = yaxis2title_ || 'Price'
        var _divid = divid_ || 'chart'
        
        var layout = {
          showlegend: true,
          legend: {orientation: 'h',
    bgcolor: 'transparent'},
          title: _title,
          yaxis: {title: _yaxistitle},
          yaxis2: {
            type: 'log',
            autorange: true,
            title: _yaxis2title,
            titlefont: {color: 'rgb(148, 103, 189)'},
            tickfont: {color: 'rgb(148, 103, 189)'},
            overlaying: 'y',
            side: 'right'
          },
          xaxis: {
            type: 'date'
          }
        };

        Plotly.newPlot(_divid, data2, layout).then(postPlotHandler('loader_c1'));
      }
 //silly little hack function to wrap a 1d array to 2d
 function wrarr(_X) {
   var X = deepClone(_X);
   var o = [];
   for(var i = 0; i<X.length;i++) {
     o[i] = [X[i]];
   }
   return(o);
 }
      
 function mOLS(_Y, _X) {
    var Y = deepClone(_Y);
    var X = deepClone(_X)
   //Multivariable linear regression
   //using linear algebra
   //first we need to add a constant to _X
   //if X is 2d matrix array [[x1,x2]...n]
   //then we need to insert another column of 1s
   // i.e. [[1, x1, x2]...n]
   var n = X.length
   var constant = Array(n).fill().map((_, i) => 1);
  // console.log(constant)
   var x = [];
   var ncol = X[0].length;
  // console.log(ncol)
   x.push(constant)
    for(var i = 0; i<ncol;i++) {
      x.push(X.map(function(value,index) { return value[i]; }));
     // console.log("adding col");
    }
    var xt = x;
    x = math.transpose(x)
    var xtx = math.multiply(xt,x); 
    var xtxinv = math.inv(xtx);
    var xtxinvxt = math.multiply(xtxinv, xt); 
    var betahat = math.multiply(xtxinvxt,Y);    
 // console.log(xtxinvxt);
    //also nice to get the hat matrix
    // (i.e. : X(XTX)−1XT)
    var hatmatrix = math.multiply(x,xtxinv)
    hatmatrix = math.multiply(hatmatrix, xt);
    var yhat = math.multiply(hatmatrix, Y);
    //leverage is the trace of the hat matrix
    var leverage = [];
    for(var h = 0; h< hatmatrix.length; h++) {
      leverage.push(hatmatrix[h][h]); //hat matrix will be square
    }
    var residuals = [];
    var SSE = 0;
    var SST = 0;
    var ybar = mean(Y.map(function(value,index) { return value[0]; }));
    for(var i = 0; i<Y.length;i++) {
        residuals.push(Y[i]-yhat[i]);
        SSE+=math.pow(residuals[i],2);
        SST+=math.pow(ybar-Y[i],2);
    }
    var MSE = SSE/n
    var cov = math.multiply(MSE,xtxinv);
    var tstats = [];
    var sebetas = [];
    for(i = 0; i<cov.length;i++) {
      sebetas.push(math.sqrt(cov[i][i]));
      tstats.push(betahat[i]/math.sqrt(cov[i][i]));
    }
    var R2 = 1-SSE/SST;
          AIC = n*math.log(SSE/n)+ (n + ncol) / (1 - (ncol + 2) / n)//Corrected AIC (Hurvich and Tsai, 1989)
          SSR = SST-SSE;
   var lm = {} 
    lm['betahat']=betahat;
    lm['yhat'] = yhat;
    lm['cov']=cov;
    lm['MSE']=MSE;
    lm['tstats']=tstats;
    lm['fitted']=yhat;//compatability
    lm['hat_matrix'] = hatmatrix;
    lm['leverage'] = leverage;
    lm['residuals'] = residuals;
    lm['se'] = sebetas;
    lm['SSE'] = SSE;
    lm['SSR'] = SSR;
    lm['AIC'] = AIC;
    lm['n'] = n;
    lm['nxvars']=ncol;
    lm['x'] = x;
    lm['y']=Y;
    lm['SST'] = SST;
    lm['R2'] = R2;
    
   return(lm);
 }
function rnd(n, p) {
    return math.format(n,p);
}

function test(lmRestricted, lmFull) {
  //test the hypothesis that whatever coefficients are left out of restricted model compared to full are jointly equal to zero
  // n is the sample size, q is the number of variables in the full model  and k is the number of restrictions imposed
  console.log("<<<<F TEST DIAGNOSTICS>>>");
  var SSER = lmRestricted.SSE;
  console.log("SSER:"+SSER)
  var SSEF = lmFull.SSE;
  console.log("SSEF:"+SSEF);
  var q = lmFull.nxvars;
  console.log("q:"+q);
  var k = q-lmRestricted.nxvars;
  console.log("k:"+k);
  var n = lmFull.n;
  console.log("n:"+n)
  console.log("<<<<END>>>>");
  //remember in mixed order integration, this has a non standard distribution. 
  return (((SSER-SSEF)/((n-k)-(n-1-q)))/(SSEF/(n-1-q)))
}

function ARDL(_Y, _X) {
 // console.log("_Y last"+math.exp(_Y[_Y.length-1]))
 var Y = deepClone(_Y)
 var X= deepClone(_X)
// X = rmNaN(X);
 //Y = rmNaN(Y);
  //basic ARDL 1,1 is
  //reg _Y _X L1._Y L2._Y
  //so need to add variables to _X assuming it is just s2f
  var y1lag = lag(Y, 1); //rewrote lag function to be more useful for ardl
  var y2lag = lag(Y, 2);

  y1lag = y1lag.slice(1,y1lag.length)
  X = X.slice(2,X.length)
  Y = Y.slice(2, Y.length)
  console.log("price last in y:" + math.exp(Y[Y.length-1]));
  var n = y2lag.length;

  var x = [];
    console.log("x length:"+X.length)
    console.log("Y1 length:" +y1lag.length)
    console.log("Y2.length:"+y2lag.length)
    console.log("Y length"+Y.length)
  x.push(X) 
  var xr = [];
//  xr.push(y1lag);
 // xr.push(y2lag);
  xr.push(X);
  
  x.push(y1lag);
  x.push(y2lag);  
  
  x = math.transpose(x)
  xr  = math.transpose(xr)
  
  m = mOLS(wrarr(Y), x);
  mred = mOLS(wrarr(Y), xr);
  console.log(mred.tstats);
  console.log(mred.cov);
  var Ftest = test(lmRestricted=mred, lmFull=m);
  
 // var lm={}
  //lm.model = m;
  //lm.residuals = 
  var mu = m.betahat[0];
  var b1 = m.betahat[1];
  var g1 = m.betahat[2];
  var g2 = m.betahat[3];
  var lr = b1[0]/(1-g1[0]-g2[0]);
  var adj = (g1[0]-1+g2[0]);
  var sr = -g1[0]; 

  //can be used to construct the STATA EC form or the non EC form:
  var model_parameters = {}
  var deltayt = mu[0]+adj*(Y[Y.length-2]-lr*X[X.length-2])-g2[0]*(Y[Y.length-2]-Y[Y.length-3])+b1[0]*(X[X.length-1]-X[X.length-2])
  console.log("mu:"+mu);
  console.log("adj:"+adj);
  console.log("Yt:" + Y[Y.length-1]);
  console.log("Yt-1:"+Y[Y.length-2]);
  console.log("Yt-2:"+Y[Y.length-3]);
  console.log("Xt:" + X[X.length-1]);
  console.log("Xt-1:"+X[X.length-2]);
  console.log("Xt-2:"+X[X.length-3]);
  
  console.log("LR:" + lr);
  console.log("Xt-1:"+X[X.length-2]);
  console.log("g2:"+g2);
  console.log("deltaYt-1:"+(Y[Y.length-2]-Y[Y.length-3]));
  console.log("b1:"+b1);
  console.log("deltaXt"+(X[X.length-1]-X[X.length-2]));
  console.log("F:"+ Ftest);
  model_parameters['F']=Ftest;
  model_parameters['mu']=mu;
  model_parameters['b1'] = b1;
  model_parameters['g1']=g1;
  model_parameters['g2']=g2;
  model_parameters['lr']=lr;  
  model_parameters['adj'] = adj; 
  model_parameters['sr'] = sr;
  model_parameters['lm'] = m;
  model_parameters['dyt'] = deltayt
  return(model_parameters);
}
 


//simple OLS (Ordinary Least Squares) function
//returns the OLS estimate for the single predictor case y~x
//depreciated. use mOLS instead.
function OLS(_Y,_X, intercept=true){
    var x = deepClone(_X);
      var y = deepClone(_Y);
  if(y.length!=x.length) throw  new Error("Y and X are of different lengths");
          var lm = {};
          var n = y.length;
          var xbar = mean(x);
          var ybar = mean(y);
          var beta_num = 0;
          var beta_den = 0;
          var beta1hat = 0;
          var beta0hat = 0;
          var yhat =new Array(y.length);
          var residuals = new Array(y.length);
          var SST = 0;
          var SSE = 0;
          var SSR = 0;
          var R2 = 0;
          var AIC = 0;
          var xx =new Array(y.length);
          var xy =new Array(y.length);
          if(!intercept) { 
            //no intercept, then
            //beta1 is just mean(x*y)/mean(x*x)
            for(var i = 0; i <y.length; i++) {
                  xy[i] = x[i]*y[i];
                  xx[i] = x[i]*x[i];
                   beta_den += Math.pow((x[i]-xbar),2);
            }
            var meanx = mean(xx);
           // console.log("mean xx:"+meanx);
            var meanxy = mean(xy);
          //  console.log("meanxy:"+meanxy);
            
            beta1hat = meanxy/meanx;
           // console.log(beta1hat);
            
          } else {
            //multiple loops required
            //:-1st to get coefficients
            //:-2nd to get residuals and fitted
            for (var i = 0; i < y.length; i++) {
              if(!Number.isNaN(x[i]) & !Number.isNaN(y[i])) { //skip nans
                beta_num += (x[i]-xbar)*(y[i]-ybar);
                beta_den += Math.pow((x[i]-xbar),2);
              }
            }
            beta1hat = beta_num/beta_den;
            beta0hat = ybar - beta1hat * xbar;
          }
                 
          for(var i =0;i<y.length;i++) {
            yhat[i]=beta0hat+beta1hat*x[i];
            residuals[i] = y[i]-yhat[i];
            SSE+=Math.pow(residuals[i],2);
            SST+=Math.pow(ybar-y[i],2);
          }

        var  sebeta1hat = Math.sqrt((SSE/(y.length-2))/beta_den);
        var  b1tstat = beta1hat/sebeta1hat

        var  sebeta0hat = Math.sqrt((SSE/(y.length-2))*(Math.pow(xbar,2)/beta_den))
        var  b0tstat = beta0hat/sebeta0hat

          AIC = n*Math.log(SSE/n)+4 //assuming k = 1 for simple case
          SSR = SST-SSE;
       //   console.log("SST:" +  SST + " SSE:"+SSE);
          R2 = 1- SSE/SST;
          lm['x'] = x;
          lm['y'] = y;
          lm['beta1hat'] = beta1hat;
          lm['beta0hat'] = beta0hat;
          lm['R2'] = R2;
          lm['SSE'] = SSE;
          lm['SST']=SST;
          lm['SSR']=SSR;
          lm['AIC']=AIC;
          lm['fitted']=yhat;
          lm['residuals']=residuals;
          lm['beta1_tstat']=b1tstat;
          lm['beta0_tstat']=b0tstat;
          //console.log(yhat);
        //  console.log("r2 "+R2);
          return lm;
}
var EGCalc = false;
var IRLCalc = false;
var EG = [];
var IR = [];
var IRL =[];
var stoIR;
var stoAG;

function getData(whichPlot,AG) {


  
  if(stoAG!=AG){
    EGCalc= false; //force recalculation of EG if changing aggregation
    EG =[];
  }
  stoAG = AG;
  
  if(stoIR!=AG) { //AG is the aggregation period dummy var
    IRLCalc = false;
    IR = [];
    IRL =[];
  }
  stoIR = AG;
//coinmetrics has excellent daily data available.
  let url  = 'btc.json';
  var stock=[];
  var price = [];
  var time=[];
  var flow=[];
  var sf = [];
  var meansf = [];
  var meanp =[];
   
   if(!isrun) {
     isrun == true
//  var output = fetch(url)
 // .then(res => res.json())
  //.then((out) => {
  if (out == null) { 
   out =$.getJSON({'url': url, 'async': false});  
  out = JSON.parse(out.responseText); 
}
    var prevStock = 0;
    var n = 0;
    for(var i in out.metricData.series) {
      time.push( out.metricData.series[i].time);
      stock.push(out.metricData.series[i].values[0]*1);
      price.push(out.metricData.series[i].values[1]*1);
      if(out.metricData.series[i].values[1] ==null) n++;
        var f = (out.metricData.series[i].values[0]-prevStock)*365.25; //daily flow
        flow.push(f) ;
        sf.push(out.metricData.series[i].values[0]/f);
      prevStock = out.metricData.series[i].values[0];
    }
  /*  for(i = 1; i<AG;i++) {
      meansf[i]=mean(sf.slice(0,i));
      meanp[i]=mean(price.slice(0,i));
    }
    for(var p = AG;p<sf.length;p++) {
      meansf[p]=mean(sf.slice(p-AG,p));
      meanp[p]=mean(price.slice(p-AG,p));
    }*/
//console.log(meansf);
    //get OLS:
    //first slice arrays to where we had prices
    time = time.slice(n,time.length);
    price = price.slice(n,price.length)//meanp.slice(n+1,meanp.length);//price.slice(n,price.length);
    sf = sf.slice(n, sf.length)//meansf.slice(n+1,meansf.length);//sf.slice(n,sf.length);
    console.log("last price actual:"+price[price.length-1]);
    var lnprice=price.map(function(e) {e = Math.log(e); return e;});
    var lnsf=sf.map(function(e) {e = Math.log(e); return e;});
var s2f = Number(sf[sf.length-1].toPrecision(2));
  lns = lnsf;
  lnp = lnprice;
  
   // const lm = OLS(lnprice, lnsf);
    //  const lm = mOLS(lnprice, wrarr(lnsf));
     if(ardl == null) ardl = ARDL(lnprice, lnsf);
     
 //  console.log(lnprice)
  //  console.log(lnsf)
  //console.log(lm)
   // var mstr = "ln(PriceUSD)="+rnd(lm.betahat[0],2)+"+"+rnd(lm.betahat[1],2)+"*ln(s2f)";
   var mstr = "<p>PriceUSD<sub>t</sub>="+rnd(math.exp(ardl.mu[0]),3)+"*s2f<sub>t</sub><sup>"+rnd(ardl.lr,3)+"</sup>";
   // mstr += "<br />R<sup>2</sup>:"+Number(lm.R2.toPrecision(2));
   // mstr += "<br />AIC:"+Number(lm.AIC.toPrecision(2));
   // mstr += "<br />t-stat coefficient:"+Number(lm.beta1_tstat.toPrecision(2));
  //  mstr += "<br />t-stat constant:"+Number(lm.beta0_tstat.toPrecision(2));
    mstr += "</p><p>S2F:"+s2f+"</p>";
    
    var ardlstr = "<p><b>ARDL form:</b> <br /> ln(Price<sub>t</sub> )= "+rnd(ardl.mu[0],3)+" + "+rnd(ardl.b1[0],3) +"ln(s2f<sub>t</sub>) + "+rnd(ardl.g1[0],3)+"ln(Price<sub>t-1</sub> )+ " + rnd(ardl.g2[0],3) + "ln(Price<sub>t-2</sub>) + &epsilon;<sub>t</sub> "
    ardlstr+="<br /> R<sup>2</sup>:"+rnd(ardl.lm.R2,4)+"</p>"
    ardlstr+="<p><b>EC form: </b><br /> &Delta;ln(Price<sub>t</sub> )="+rnd(ardl.mu[0],3)+" + "+rnd(ardl.adj,3)+"*[ln(Price<sub>t-1</sub>) - "+rnd(ardl.lr,3)+"*ln(s2f<sub>t-1</sub>)] - "+rnd(ardl.g2[0],3)+"*&Delta;ln(Price<sub>t-1</sub> )+ "+rnd(ardl.b1[0],3)+"*&Delta;ln(s2f<sub>t</sub>) + &epsilon;<sub>t</sub></p>"
    ardlstr+="<p>&Delta;ln(Price) = "+rnd(ardl.dyt,4) + "<br />"
  //  ardlstr+="&rarr; &Delta;Price = "+rnd(math.exp(ardl.dyt),4)+"</p>"
//console.log(ardlstr);

    document.getElementById('model').innerHTML=mstr;
    document.getElementById('ardlmodel').innerHTML=ardlstr;


   // var residuals = ardl.lm.residuals
//ok, now we need to find the likelihood of the current
//residual (eg residual[residual.length]) given
//the previous residuals. The likelihood is found by the

 meanres = mean(ardl.lm.residuals);
 stdev = jStat.stdev(ardl.lm.residuals,true);
 p05 = percentile(ardl.lm.residuals,0.05);//jStat.normal.inv( 0.05, meanres, stdev)
 p95 = percentile(ardl.lm.residuals,0.95);//jStat.normal.inv( 0.95, meanres, stdev)

 p05_L = lowerCI(0.05,0.05,ardl.lm.residuals);
 p95_U = upperCI(0.05, 0.95, ardl.lm.residuals);






  for(var i = 0; i< ardl.lm.residuals.length; i++){
       fifth.y[i]=p05;
      q05LCI.y[i] = p05_L;
      q95UCI.y[i] = p95_U;
      q05LCIexp.y[i] = Math.exp(ardl.lm.yhat[i]+p05_L);
      q95UCIexp.y[i] = Math.exp(ardl.lm.yhat[i]+p95_U)
      mean_.y[i]=meanres;
      ninetyfifth.y[i]=p95;

      q95.y[i]=Math.exp(ardl.lm.yhat[i]+p95);
      q05.y[i]=Math.exp(ardl.lm.yhat[i]+p05);
      fitted.y[i]=ardl.lm.yhat[i];
      rvf.x[i]=ardl.lm.yhat[i];
      rvf.y[i]= ardl.lm.residuals[i];
    }
    likelihood.y = getLikelihood(ardl.lm.residuals,0);
    console.log("likelihyood leng:"+likelihood.y.length);
    
    likelihood.x = time.map(function(value,index) {
      var d = new Date(value);
      
      var datestring = + d.getFullYear() +"-" + ("0"+(d.getMonth()+1)).slice(-2) +"-"+ ("0" + d.getDate()).slice(-2)
      return datestring });
      likelihood.x = likelihood.x.slice(2,likelihood.x.length);
   console.log(" date length:"+likelihood.x.length)
    residuals.x = likelihood.x;
          residuals.y=ardl.lm.residuals;
       //fitted.y= lm.yhat.map(function(e) {e = Math.exp(e); return e;});
var prices={y: price.slice(2, price.length), x:likelihood.x, name: 'price', type:'scatter', yaxis:'y2', mode:'markers', marker: {color:likelihood.y}};
//var prettydata = [prettyprices] 
          fitted.x=residuals.x;
    //prettyprices.x=residuals.x;
    //prices.y=price;
    q95.x = prices.x;
    q05.x = prices.x;
    q05LCI.x = prices.x;
    q95UCI.x = prices.x;
    q95UCIexp.x = prices.x;
    q05LCIexp.x = prices.x;
    fifth.x = prices.x;
    mean_.x = prices.x;
    ninetyfifth.x = prices.x;
  
     ll = likelihood.y[likelihood.y.length-1];
 //console.log(ll);
  lr =    ardl.lm.yhat[ardl.lm.yhat.length-1];
 lr = lr[0];
  pr = prices.y[prices.y.length-1];
  diff = pr - math.exp(lr);
  meandistance = (residuals.y[residuals.y.length-1]);
  tomorrow = math.exp(ardl.mu[0]+ardl.g1[0]*lnprice[lnprice.length-1]+ardl.g2[0]*lnprice[lnprice.length-2]+ardl.b1[0]*lnsf[lnsf.length-1])
var valued = "fair"
 if(diff >0) {
  if(ll<0.25) valued = "fair to very mildly overvalued";
  if(ll<0.15) valued = "mildly overvalued";
  if(ll<0.10) valued = "overvalued";
  if(ll<0.05) valued = "extremely overvalued";
 } else {
  if(ll<0.25) valued = "fair to very mildly undervalued";
  if(ll<0.15) valued = "mildly undervalued";
  if(ll<0.10) valued = "undervalued";
  if(ll<0.05) valued = "extremely undervalued";
}
var ovund = ""
var kelly = ((0.5-ll)+0.5)*(1+1)-1
if(diff <0) ovund = "under"
if(diff >=0) ovund = "over"
var tabled = "Todays Bitcoin price of $"+rnd(pr,4)+" is  $" + rnd(diff,2) + " "+ ovund+" the ARDL estimation of $"+rnd(math.exp(lr),4)+". The likelihood of this occuring is about "+rnd(ll,2)+". "+
"Therefore according to the ARDL stock to flow model, the Bitcoin price is " + valued+". <p>The Kelly fraction assuming equal odds is "+rnd(kelly,2) + ".</p>"
//tabled+="<p>According to the ARDL; given todays price, tomorrows should be around:$"+rnd(tomorrow,4)+"</p>"
if(document.getElementById("prettychart")!=null) document.getElementById("prettychart").innerHTML=tabled;

var layout = {
  width: 500,
  height: 400,
  margin: { t: 25, r: 25, l: 25, b: 25 },
  paper_bgcolor: "transparent",
  font: { color: "darkblue", family: "Arial" }
};


var datarl = [
  {
    type: "indicator",
    mode: "gauge+number",
    value: ll,
    title: { text: "Residual Likelihood", font: { size: 24 } },
    gauge: {
      axis: { range: [null, 1], tickwidth: 1, tickcolor: "darkblue" },
      bar: { color: "darkblue" },
      bgcolor: "white",
      borderwidth: 2,
      bordercolor: "gray",
      steps: [
        { range: [0, 0.05], color: "red" },
        { range: [0.05, 0.15], color: "orange"},
        { range: [0.15, 1], color: "green"}
      ]
    }
  }
];

var layout = {
  width: 500,
  height: 400,
  margin: { t: 25, r: 25, l: 25, b: 25 },
  paper_bgcolor: "transparent",
  font: { color: "darkblue", family: "Arial" }
};



var datadm = [
  {
    type: "indicator",
    mode: "gauge+number+delta",
    value: rescale(1.96,-1.96,p95_U, p05_L,meandistance),
    title: { text: "Distance from Mean", font: { size: 24 } },
    delta: { reference: 0 },
    gauge: {
      axis: { range: [rescale(1.96,-1.96,p95_U, p05_L,p05_L*2),rescale(1.96,-1.96,p95_U, p05_L,p95_U*2)], tickwidth: 1, tickcolor: "darkblue" },
      bar: { color: "darkblue" },
      bgcolor: "white",
      borderwidth: 2,
      bordercolor: "gray",
      steps: [
        { range: [rescale(1.96,-1.96,p95_U, p05_L,p05_L*2), rescale(1.96,-1.96,p95_U, p05_L,p05_L)], color: "green" },
        { range: [rescale(1.96,-1.96,p95_U, p05_L,p95_U), rescale(1.96,-1.96,p95_U, p05_L,p95_U*2)], color: "red"}
      ],
      threshold: {
        line: { color: "red", width: 4 },
        thickness: 0.75,
        value: 0
      }
    }
  }
];



  //});


  }
    if(whichPlot==1) {
      var residualz={y: likelihood.y, x:residuals.x, name: 'residualz', type: 'scatter', yaxis: 'y1',  mode: 'line', marker:{color:residuals.y}};

var datax = [residualz];

     csvContent="date,price, residuals, likelihood\n";
       for(j =0;j<prices.y.length;j++) {
         csvContent+=prices.x[j]+","+ prices.y[j]+","+residuals.y[j]+","+ likelihood.y[j]+"\n";
       }
plot1(datax,"Likelihood");

} else if (whichPlot==2) {
  var data3 = [prices,q05LCI, fifth, ninetyfifth, q95UCI, residuals, mean_];
  

     csvContent="date,price,q05LCI, fifth, ninetyfifth, q95UCI, residuals, mean\n";
       for(j =0;j<prices.y.length;j++) {
         csvContent+=prices.x[j]+","+ prices.y[j]+","+q05LCI.y[j]+","+ fifth.y[j]+","+  ninetyfifth.y[j]+","+ q95UCI.y[j]+","+ residuals.y[j]+","+  mean_.y[j]+"\n";
       }

 
plot2(data3)
} else if (whichPlot==3) {
 
 
 //console.log("price:"+pr+" exonentiatiated model:"+Math.exp(Math.log(pr)-lr)+" diff:"+diff)

 Plotly.newPlot('gg', datarl, layout)//loader_c1 clears when div loads anyway .then(postPlotHandler('loader_c1'));

Plotly.newPlot('gr', datadm, layout);


//dfullert=egranger(lnprice, lnsf);
/*
var dataci = [
  {
    type: "indicator",
    mode: "gauge+number",
    value: dfullert,
    title: { text: "Cointegration", font: { size: 24 } },
    gauge: {
      axis: { range: [-24,0], tickwidth: 1, tickcolor: "darkblue" },
      bar: { color: "darkblue" },
      bgcolor: "white",
      borderwidth: 2,
      bordercolor: "gray",
      steps: [
        { range: [-24, -3.5], color: "green" },
        { range: [-3.5, -2.5], color: "orange"},
        { range: [-2.5, 0], color:"red"}
      ],

    }
  }
];



Plotly.newPlot('coint_gauge', dataci, layout);
} else if (whichPlot==4) {
    var data = [q05LCIexp, q05, q95,q95UCIexp,prices];
    csvContent="date,price, q05LCI, q05, q95,q95UCI\n";
    for(j =0;j<prices.y.length;j++) {
         csvContent+=prices.x[j]+","+ prices.y[j]+","+q05LCIexp.y[j]+","+ q05.y[j]+","+ q95.y[j]+","+ q95UCIexp.y[j]+"\n";
    }
    var layout = {
         title: 'S2F BTC Model',
         yaxis2: {
           type: 'log',},
         xaxis: {
           type: 'date'
         }, 
         showlegend: true,
        legend: {orientation: 'h',
        bgcolor: 'white'}
       };
    Plotly.newPlot('model_chart', data, layout).then(postPlotHandler('loader_c1'));;
   
   if (EGCalc == false) {
    var r = confirm("Warning: the Engle-Granger stat over time will take a long time to compute. Continue?")
    if(r ==true) {
      for(i = 2; i<time.length;i++){
        EG[i]= egranger(lnprice.slice(0, i), lnsf.slice(0,i));
      } 
      EGCalc = true;
    }
  } 
    var eg_data=[{y: EG, x:time, name: 'Engle-Granger', type: 'scatter', yaxis: 'y1',  mode: 'lines'}];


  var layouteg = {
      title:'Engle-Granger Statistic',
      showLegend: true, 
      legend: {orientation:'h',
      bgcolor: 'white'}};
      Plotly.newPlot('eg_chart', eg_data, layouteg).then(postPlotHandler('loader_c1'));

    csvContentEG="date, EngleGranger\n";
       for(j =0;j<time.length;j++) {
         csvContentEG+=time[j]+","+ EG[j]+"\n";
       }
    */
} else if (whichPlot == 5) {   
  var data2 = [rvfPlot(ardl.lm)];

  //console.log(data2)
  var layout2 = {
         title: 'Residuals v Fitted', 
         showlegend: true,
        legend: {orientation: 'h',
    bgcolor: 'white'}
       };
    Plotly.newPlot('rvf_chart', data2, layout2).then(postPlotHandler('loader_c1'));;
    
    var data3 = [hist(ardl.lm)];
    var layoutHis = {
      title: 'Histogram of Residuals',
      showlegend: false,
      bgolor:'white'};
      Plotly.newPlot('hist_chart', data3, layoutHis).then(postPlotHandler('loader_c12'));
  var data4 = [linearity(lnsf, lnprice)];
  var layoutLin = {
      title:'Linearity',
      showLegend: true, 
      legend: {orientation:'h',
      bgcolor: 'white'}};
      Plotly.newPlot('lin_chart', data4, layoutLin).then(postPlotHandler('loader_c13'));

} else if (whichPlot == 6) {

  
   if (IRLCalc == false) {
    var res = confirm("Warning: the Instantaneous Residual Likelihood over time will take a long time to compute. Continue?")
    if(res ==true) {
      for(i = 30; i<time.length;i++){
        //first, get OLS for each point i
        var am = ARDL(lnprice.slice(0,i), lnsf.slice(0,i));
        var m = am.lm;
        //now save residuals for point i
        IR[i]=m.residuals[i-1];
       // console.log(IR[i-1]);
      }
      //now get likelhoods
      IRL = [NaN,NaN].concat(getLikelihood(IR,0));
     // console.log(IRL);
      IRLCalc = true;
    }
  } 
  var irlprices={y: prices.y, x:prices.x, name: 'price', type:'scatter', yaxis:'y2', mode:'markers',  marker: {color:IR}};
    var irl_plot={y: IRL, x:prices.x, name: 'IRL', type: 'scatter', yaxis: 'y1',  mode: 'lines'};
    var irldata = [irl_plot, irlprices];
 var irllayout = {
          showlegend: true,
          legend: {orientation: 'h',
    bgcolor: 'transparent'},
          title: 'Instantaneous Residual Likelihood Analysis',
          yaxis: {title: 'Residual Likelihood'},
          yaxis2: {
            type: 'log',
            autorange: true,
            title: 'Price USD',
            titlefont: {color: 'rgb(148, 103, 189)'},
            tickfont: {color: 'rgb(148, 103, 189)'},
            overlaying: 'y',
            side: 'right'
          },
          xaxis: {
            type: 'date'
          }
        };

        Plotly.newPlot('irlchart', irldata, irllayout).then(postPlotHandler('loader_c1'));
  //  plot2(irldata,"Immediate Residuals", "Residuals", "Price","irlchart")
  /*  
  var layoutirl = {
      title:'Immediate Residual Likelihood',
      showLegend: true, 
      legend: {orientation:'h',
      bgcolor: 'white'}};
      Plotly.newPlot('irlchart', irl_plot, layoutirl).then(postPlotHandler('loader_c1'));
*/
    csvContentIRL="date, IRL\n";
       for(j =0;j<time.length;j++) {
         csvContentIRL+=time[j]+","+ IRL[j]+"\n";
       }
    
}  

}
 

function activateLink(linked) {
  $('.nav').each(function(i, obj) {
    obj.css=({"float": "left;", "display": "block;","color": "#f2f2f2;","text-align": "center;","padding": "14px 16px;","text-decoration":"none;"});
});
$(linked).css=({"background-color": "#ddd;","color": "black;"})

}



function rvfPlot( _lm_model) {
    var lm_model = JSON.parse(JSON.stringify(_lm_model))
    var fitted = [];
    var res = lm_model.residuals;
    var y = lm_model.y;
    for(var i = 0; i<res.length; i++) {
        fitted[i] = y[i]+res[i]
    }
var rvf={y: res.slice(1,res.length), x:fitted.slice(1,fitted.length), name: 'RVF', type: 'scatter', yaxis: 'y1',  mode: 'markers'};

   return rvf;
}

function hist( _lm_model) {
    var lm_model = JSON.parse(JSON.stringify(_lm_model))
  var res = lm_model.residuals;
  var trace = {
      x: res,
      type: 'histogram',
    };
    return trace;
}

function linearity(lns2f,lnprice) {
  //var lm_model = JSON.parse(JSON.stringify(_lm_model))
   // var price = lm_model.y[0];
  //  var s2f = math.transpose(lm_model.x)[1];
    var trace = {y: lnprice, x: lns2f, name: "Linearity", type: 'scatter', mode:'markers', yaxis:'y1'};
    return trace;
}


// The download function takes a CSV string, the filename and mimeType as parameters
// Scroll/look down at the bottom of this snippet to see how download is called
function downloadz(content, fileName, mimeType) {
  var a = document.createElement('a');
  mimeType = mimeType || 'application/octet-stream';

  if (navigator.msSaveBlob) { // IE10
    navigator.msSaveBlob(new Blob([content], {
      type: mimeType
    }), fileName);
  } else if (URL && 'download' in a) { //html5 A[download]
    a.href = URL.createObjectURL(new Blob([content], {
      type: mimeType
    }));
    a.setAttribute('download', fileName);
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  } else {
    location.href = 'data:application/octet-stream,' + encodeURIComponent(content); // only this mime type is supported
  }
}



