Simulation = require('./cacatoo.js') // Loads the Simulation class from installation, or from local package like below

let yargs = require('yargs')
let cmd_params = yargs.argv

// Import paramter values from command line
// to run on the command line: node model.js --pi 4 --R 1 --seed 1 --VD -1
let pi = 10**(cmd_params.pi * -0.5)
let radius = cmd_params.R   // If =1, R is set to infinity
let iter = cmd_params.seed
let VD = cmd_params.VD

let Vdiff = 10**(-1*VD) // If >1 (i.e. VD < 0), the virions are well-mixed
let seed = iter

// Parameters

b = 0.5
phi = 1             // interaction strength parameter (also carrying capacity parameter)
d = b*phi/(1+phi)   // max increase in death rate
alpha = 0.0005       // lysis rate
dint = 0.005

mu = 0.0005       // gene loss rate
loss = 0.0005    // phage loss rate

if (pi < 10**-5) pi = 0.0 // Asume no privatization if pi is below 0.00001

beta = 0.005  // infection rate
rV = 30      // burst size
g = 0.005     // virion decay rate

var vir = ["Uc", "Lc", "Lp", "Lcp"] // Define list of carrier genotypes

var size = 120
let num_nT = 0, num_T = 0

Number.prototype.mod = function(n) {
    return ((this%n)+n)%n;
}
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Configuration
let config = {
    title: "Within host bacteria-phage dynamics",
    description: "Radius = 0, phi = 1",
    maxtime: 250000,
    seed: seed,
    ncol: size,
    nrow: size,		            // dimensions of the grid to build
    wrap: [true, true]       // Wrap boundary [COLS, ROWS]
}

sim = new Simulation(config)

// Creating two overlaying grids - one for the local environments with virions and the other for the bacteria
sim.makeGridmodel("cell");
sim.makeGridmodel("ext");

// Initialization of bacteria with all Lc genotype
sim.cell.initialise = function () {
    sim.initialGrid(sim.cell, "alive", 0, 0.0, 'U', 0, 'Uc', 0, 'L', 0.0, 'Lc', 1.0, 'Lp', 0.0, 'Lcp', 0.0)
}

sim.cell.initialise()

// Initialization of external virion particle density
sim.ext.initialise = function () {
    sim.initialGrid(sim.ext, "VT", 0, 1)
    sim.initialGrid(sim.ext, "V", 0, 1)
}

sim.ext.initialise()

// Next State Function of bacteria
sim.cell.nextState = function (i, j) {

    var V = sim.ext.grid[i][j].V
    var VT = sim.ext.grid[i][j].VT

    // If gp is empty, copy a random neighbour with probability b
    if (this.grid[i][j].alive == 0){
        var rand = sim.rng.genrand_real1()
        if (rand < b){
            nbr = this.randomMoore8(this, i, j)
            this.grid[i][j].alive = nbr.alive
        }
    }

    // If gp is occupied, cell dies or lyses (by induction). If it survives, it may get infected or mutate. Only one of these events happen
    else {

        // If R>1, the fraction of virulent strains in the circular nbd is found
        if (radius>1){
          var numL_nT = 0 // local number of non-carrier cells
          var numL_T = 0 // local number of carrier cells

          for (let x = i-radius; x < i+radius; x++){                          // x are columns
              for (let y = j-radius; y < j+radius; y++){                           // y are rows
                  if ((Math.pow((i - x), 2) + Math.pow((j - y), 2)) < Math.pow(radius, 2)){
                      if (vir.includes(this.grid[x.mod(this.nc)][y.mod(this.nr)].alive)) numL_T ++
                      else if (this.grid[x.mod(this.nc)][y.mod(this.nr)].alive != 0) numL_nT ++
                  }
              }
           }
        }

        // If R=1, local numbers are set to global numbers (i.e. radius is infinite)
        else{
            var numL_nT = num_nT
            var numL_T = num_T
        }

        // If cell is U, die or get infected or do nothing
        if (this.grid[i][j].alive == 'U'){
            // death
            frac = numL_T/(numL_nT+numL_T)
            var rand = sim.rng.genrand_real1()
            var death = dint + d*(1 - (1-pi)*frac)
            if (rand < death){
                this.grid[i][j].alive = 0
            }
            // or get infected
            else if (rand < death + beta*(V+VT)){
                // if infected, choose virion genotype
                var rand1 = sim.rng.genrand_real1()
                if (rand1 < V/(V+VT)){
                    this.grid[i][j].alive = 'L'
                    sim.ext.grid[i][j].V -= 1
                }
                else {
                    this.grid[i][j].alive = 'Lp'
                    sim.ext.grid[i][j].VT -= 1
                }
            }

        }

        // If cell is Uc, die or get infected or lose gene or do nothing
        if (this.grid[i][j].alive == 'Uc'){
            // death
            frac = (numL_T)/(numL_nT+numL_T)
            var rand = sim.rng.genrand_real1()
            var death = dint + d*(1 - pi - (1-pi)*frac)
            if (rand < death){
                this.grid[i][j].alive = 0
            }

            else if (rand < death +  mu){
                    this.grid[i][j].alive = 'U'
            }
            else if (rand < death + mu + beta*(V+VT)){
                var rand1 = sim.rng.genrand_real1()
                if (rand1 < V/(V+VT)){
                    this.grid[i][j].alive = 'Lc'
                    sim.ext.grid[i][j].V -= 1
                  }
                else {
                    this.grid[i][j].alive = 'Lcp'
                    sim.ext.grid[i][j].VT -= 1
                }
            }

        }

        // If cell is L, die or lyse by induction or lose phage or do nothing
        if (this.grid[i][j].alive == 'L'){
            // death
            frac = numL_T/(numL_nT+numL_T)
            var rand = sim.rng.genrand_real1()
            var death = dint + d*(1 - (1-pi)*frac)
            if (rand < death){
                this.grid[i][j].alive = 0
            }
            else if(rand  < death + alpha){
                    this.grid[i][j].alive = 0
                    sim.ext.grid[i][j].V += rV
            }
            else if (rand < death + alpha + loss){
                    this.grid[i][j].alive = 'U'
            }
        }

        // If cell is Lc, die or lyse by induction or lose phage or lose gene or do nothing
        if (this.grid[i][j].alive == 'Lc'){
            // death
            frac = (numL_T)/(numL_nT+numL_T)
            var rand = sim.rng.genrand_real1()
            var death = dint + d*(1 - pi -(1-pi)*frac)
            if (rand < death){
                this.grid[i][j].alive = 0
            }
            else if(rand  < death + alpha){
                    this.grid[i][j].alive = 0
                    sim.ext.grid[i][j].V += rV
            }
            else if (rand < death + alpha + loss){
                    this.grid[i][j].alive = 'Uc'
            }
            else if (rand < death + alpha + loss + mu){
                    this.grid[i][j].alive = 'L'
            }
        }

        // If cell is Lp, die or lyse by induction or lose phage or lose gene or do nothing
        if (this.grid[i][j].alive == 'Lp'){
            // death
            frac = (numL_T)/(numL_nT+numL_T)
            var rand = sim.rng.genrand_real1()
            var death = dint + d*(1 - pi -(1-pi)*frac)
            if (rand < death){
                this.grid[i][j].alive = 0
            }
            else if(rand  < death + alpha){
                    this.grid[i][j].alive = 0
                    sim.ext.grid[i][j].VT += rV
            }
            else if (rand < death + alpha + loss){
                    this.grid[i][j].alive = 'U'
            }
            else if (rand < death + alpha + loss + mu){
                    this.grid[i][j].alive = 'L'
            }
        }

        // If cell is Lcp, die or lyse by induction or lose phage or lose gene on either locus or do nothing
        if (this.grid[i][j].alive == 'Lcp'){
            // death
            frac = (numL_T)/(numL_nT+numL_T)
            var rand = sim.rng.genrand_real1()
            var death = dint + d*(1 - pi -(1-pi)*frac)
            if (rand < death){
                this.grid[i][j].alive = 0
            }
            else if(rand  < death + alpha){
                    this.grid[i][j].alive = 0
                    sim.ext.grid[i][j].VT += rV
            }
            else if (rand < death + alpha + loss){
                    this.grid[i][j].alive = 'Uc'
            }
            else if (rand < death + alpha + loss + mu){
                    this.grid[i][j].alive = 'Lc'
            }
            else if (rand < death + alpha + loss + 2*mu){
                    this.grid[i][j].alive = 'Lp'
            }
        }

    }

}

// Next State Function of virion particle density
sim.ext.nextState = function (i, j) {
    // if Vdiff <= 1 virions can diffuse
    if (Vdiff <= 1){
      // apply diffusion and decay with probability to each non-carrier virion particle
      var moved = 0
      var decay = 0
      for (let k = 0; k < this.grid[i][j].V; k++) {
          if (sim.rng.genrand_real1() < g){
              decay ++
          }
          else if (sim.rng.genrand_real1() < Vdiff){
              nbr = this.randomMoore8(this, i, j)
              nbr.V ++
              moved ++
          }
      }
      this.grid[i][j].V -= moved + decay

      // apply diffusion and decay with probability to each carrier virion particle
      var moved = 0
      var decay = 0
      for (let k = 0; k < this.grid[i][j].VT; k++) {
          if (sim.rng.genrand_real1() < g){
              decay ++
          }
          else if (sim.rng.genrand_real1() < Vdiff){
              nbr = this.randomMoore8(this, i, j)
              nbr.VT ++
              moved ++
          }
      }
      this.grid[i][j].VT -= moved + decay
    }

    // If Vdiff >1, virions are well-mixed
    else {
      var moved = 0
      var decay = 0
      for (let k = 0; k < this.grid[i][j].V; k++) {
          if (sim.rng.genrand_real1() < g){
              decay ++
          }
          else {
              var p = sim.rng.genrand_int(0,this.nc-1)
              var q = sim.rng.genrand_int(0,this.nr-1)
              nbr = this.grid[p][q]
              nbr.V ++
              moved ++
          }
      }
      this.grid[i][j].V -= moved + decay

      var moved = 0
      var decay = 0
      for (let k = 0; k < this.grid[i][j].VT; k++) {
          if (sim.rng.genrand_real1() < g){
              decay ++
          }
          else {
              var p = sim.rng.genrand_int(0,this.nc-1)
              var q = sim.rng.genrand_int(0,this.nr-1)
              nbr = this.grid[p][q]
              nbr.VT ++
              moved ++
          }
      }
      this.grid[i][j].VT -= moved + decay
    }
}


sim.cell.update = function (){

    // Compute abundances by cell-type
    num_nT = 0, num_T = 0
    let Utot = 0, Uctot = 0, Ltot = 0, Lctot = 0, Lptot = 0, Lcptot = 0, Vtot = 0, VTtot = 0
    for (let i = 0; i < this.nc; i++) for (let j = 0; j < this.nr; j++) {
        if (vir.includes(this.grid[i][j].alive)) num_T ++
        else if (this.grid[i][j].alive != 0) num_nT ++

        if (this.grid[i][j].alive == 'U') Utot ++
        else if (this.grid[i][j].alive == 'Uc') Uctot ++
        else if (this.grid[i][j].alive == 'L') Ltot ++
        else if (this.grid[i][j].alive == 'Lc') Lctot ++
        else if (this.grid[i][j].alive == 'Lp') Lptot ++
        else if (this.grid[i][j].alive == 'Lcp') Lcptot ++

        Vtot += sim.ext.grid[i][j].V
        VTtot += sim.ext.grid[i][j].VT
    }

    this.synchronous()
    //this.perfectMix()

    // Write output
    if (sim.time%1000==0){
      sim.write_append(sim.time+'\t'+Utot+'\t'+Uctot+'\t'+Ltot+'\t'+Lctot+'\t'+Lptot+'\t'+Lcptot+'\t'+Vtot+'\t'+VTtot+'\n', "timeseries_data/sample/sample_R"+cmd_params.R+"pi"+cmd_params.pi+"iter"+iter+"VD"+VD+".dat")
      if (Uctot+Lctot+Lptot+Lcptot+VTtot==0) process.exit() // end simulation is gene is extinct
    }

    // induce mutation (1% Lcbecomes Lp) is phage encoded gene goes extinct
    if (Lptot + Lcptot + VTtot == 0){
      for (let i = 0; i < this.nc; i++) for (let j = 0; j < this.nr; j++) {
        if (this.grid[i][j].alive == 'Lc' && sim.rng.random() < 0.01){
          this.grid[i][j].alive = 'Lp'
        }
      }
    }
}

sim.ext.update = function (){
    this.asynchronous()
}

sim.start()
