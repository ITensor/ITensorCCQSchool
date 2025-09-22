using ITensors
using ITensorMPS

using ITensors, ITensorMPS

let
  #Build the physical indices for 30 spins (spin 1/2)
  N = 30
  sites = siteinds("S=1/2",N)

  #Build the Heisenberg Hamiltonian as an MPO
  os = OpSum()
  for j=1:N-1
    os += "Sz",j,"Sz",j+1
    os += 1/2,"S+",j,"S-",j+1
    os += 1/2,"S-",j,"S+",j+1
  end
  H = MPO(os,sites)

  #It has bond dimension 5
  println("MPO bond dimension is $(ITensorMPS.maxlinkdim(H))") 

  #Initial state for DMRG
  psi0 = random_mps(sites;linkdims=10)

  #DMRG Parameters
  nsweeps = 5
  maxdim = [10,20,100,100,200]
  cutoff = [1E-10]

  #Run DMRG!
  energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)

  return
end