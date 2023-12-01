boltz(V::Float64,A::Float64,B::Float64) = 1/(1 + exp((V+A)/B))
tauX(V::Float64,A::Float64,B::Float64,D::Float64,E::Float64) = A - B/(1+exp((V+D)/E))
mNainf(V::Float64) = boltz(V,35.5,-5.29)
taumNa(V::Float64) = tauX(V,1.32,1.26,120.,-25.)
hNainf(V::Float64) = boltz(V,48.9,5.18)
tauhNa(V::Float64) = (0.67/(1+exp((V+62.9)/-10.0)))*(1.5 + 1/(1+exp((V+34.9)/3.6)))
mKdinf(V::Float64) = boltz(V,12.3,-11.8)
taumKd(V::Float64) = tauX(V,7.2,6.4,28.3,-19.2)
mCaTinf(V::Float64) = boltz(V,67.1,-7.2)
taumCaT(V::Float64) = tauX(V,21.7,21.3,68.1,-20.5)
hCaTinf(V::Float64) = boltz(V,80.1,5.5)
tauhCaT(V::Float64) = 2*tauX(V,205.,89.8,55.,-16.9)
mHinf(V::Float64) = boltz(V,80.,6.)
taumH(V::Float64) = tauX(V,272.,-1149.,42.2,-8.73)
mKCainf(Ca::Float64) = (Ca/(Ca+Kd))^2
Tm(V::Float64) = 1/(1+exp(-(V-2)/5))

ICaT(V::Float64,mCaT::Float64,hCaT::Float64,gCaT::Float64) = gCaT*mCaT^3*hCaT*(V-VCa)

function dV(gNa::Float64, gKd::Float64, gl::Float64, V::Float64, mNa::Float64, hNa::Float64, mKd::Float64, mCaT::Float64, hCaT::Float64, mH::Float64, Ca::Float64, Iapp::Float64, Istep::Float64, gCaT::Float64, gH::Float64, gKCa::Float64)
  (dt)*(1/C)*(-gNa*mNa^3*hNa*(V-VNa) -gKd*mKd^4*(V-VK) -gl*(V-Vl) -ICaT(V,mCaT,hCaT,gCaT) -gKCa*mKCainf(Ca)*(V-VK) -gH*mH*(V-VH) +Iapp +Istep)
end

dmNa(V::Float64,mNa::Float64) = (dt)*((1/taumNa(V))*(mNainf(V) - mNa))
dhNa(V::Float64,hNa::Float64) = (dt)*((1/tauhNa(V))*(hNainf(V) - hNa))
dmKd(V::Float64,mKd::Float64) = (dt)*((1/taumKd(V))*(mKdinf(V) - mKd))
dmCaT(V::Float64,mCaT::Float64) = (dt)*((1/taumCaT(V))*(mCaTinf(V) - mCaT))
dhCaT(V::Float64,hCaT::Float64) = (dt)*((1/tauhCaT(V))*(hCaTinf(V) - hCaT))
dmH(V::Float64,mH::Float64) = (dt)*((1/taumH(V))*(mHinf(V) - mH))
dCa(V::Float64,mCaT::Float64,hCaT::Float64,Ca::Float64,gCaT::Float64,k1::Float64,k2::Float64) = (dt)*(-k1*ICaT(V,mCaT,hCaT,gCaT) -k2*Ca)
dAMPA(V::Float64,AMPA::Float64) = (dt)*(1.1*Tm(V)*(1-AMPA)-0.19*AMPA)
dGABAA(V::Float64,GABAA::Float64) = (dt)*(0.53*Tm(V)*(1-GABAA)-0.18*GABAA)
dGABAB(V::Float64,GABAB::Float64) = (dt)*(0.016*Tm(V)*(1-GABAB)-0.0047*GABAB)


function dx(x::Float64, V::Float64, Vprev::Float64)
  (dt)*-x/tau_p + 1*(Vprev<= 0 && V >0)
end

function dx2(x::Float64, V::Float64, Vprev::Float64)
  (dt)*-x/tau_x + 1*(Vprev<= 0 && V >0)
end

function dy(y::Float64, V::Float64, Vprev::Float64)
  (dt)*-y/tau_m + 1*(Vprev<= 0 && V >0)
end

function dy2(y::Float64, V::Float64, Vprev::Float64)
  (dt)*-y/tau_y + 1*(Vprev<= 0 && V >0)
end


function dw(new::Float64, w::Float64, tau::Float64)
  (dt)*(new - w)*1/tau
end

function dl(w_dot::Float64, tau_l::Float64) 
  (dt)*(-w_dot)*1/tau_l
end 

function dl_tonic(w_dot::Float64, tau_l::Float64)
  (dt)*(w_dot)*1/tau_l
end 




function simulate_scenario_noisy(freq::Float64, delay_pre::Float64, w::Float64, l::Float64)
  #SYN = [(gEE/4)*(rand(nEcells,nEcells)-0.5*ones(nEcells, nEcells)).+gEE (gEI/4)*(rand(nEcells,nIcells)-0.5*ones(nEcells, nIcells)).+gEI;(gIE/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE (gII/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII;(gIE2/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE2 (gII2/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII2]
  gEC = gEC_const*w

  # Synaptic connections
  #SYN = [(gEE/4)*(rand(nEcells,nEcells)-0.5*ones(nEcells, nEcells)).+gEE (gEI/4)*(rand(nEcells,nIcells)-0.5*ones(nEcells, nIcells)).+gEI;(gIE/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE (gII/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII;(gIE2/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE2 (gII2/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII2]
  SYN = [gEE*ones(nEcells,nEcells) gEI*ones(nEcells,nIcells) gEC*ones(nEcells, nCcells);gIEGABAA*ones(nIcells,nEcells) gIIGABAA*ones(nIcells,nIcells) gICGABAA*ones(nIcells, nCcells);gIEGABAB*ones(nIcells,nEcells) gIIGABAB*ones(nIcells,nIcells) gICGABAB*ones(nIcells, nCcells);gCE*ones(nCcells,nEcells) gCI*ones(nCcells,nIcells) gCC*ones(nCcells,nCcells)]
  print(SYN[1,3],' ' , SYN[2,3])
  # Initial conditions
  V=-60*ones(ncells)
  Vprev=-60*ones(ncells)
  x::Float64 = 0.
  y::Float64 = 0.
  y2::Float64 = 0.
  x2::Float64 = 0.

  x_prev::Float64 = 0.
  y_prev::Float64 = 0.
  stop = 0
  # VV = zeros(Tdt, ncells)
  VV=zeros(convert(Int64,Tdt*dt), ncells)
  ww=zeros(convert(Int64,Tdt*dt))
  ll=zeros(convert(Int64,Tdt*dt))


  idx_w=1

  mNa=mNainf(V[1])*ones(ncells)
  hNa=hNainf(V[1])*ones(ncells)
  mKd=mKdinf(V[1])*ones(ncells)
  mCaT=mCaTinf(V[1])*ones(ncells)
  hCaT=hCaTinf(V[1])*ones(ncells)
  mH=mHinf(V[1])*ones(ncells)
  Ca=((-k1_E/k2_E)*ICaT(V[1],mCaT[1],hCaT[1],gCaT_E[1]))*ones(ncells)
  AMPA=zeros(ncells)
  GABAA = zeros(ncells)
  GABAB = zeros(ncells)
  last = 2
  steady = 0.
  A = 0.
  new = w
  tau = 0.1
  count_pre = 0
  count_post = 0
  w_dot=0


  

  zdt_freqE = 0
  z_delayed = delay_pre
  X = Distributions.Gaussian(freq,0.3*freq)
  freqE = rand(X)
  spikedE = 0
  spikedC = 0
  newFreqFlag = 0
  newfreqE = 0

  burst::Int64 = 0
  t_Spike = []
  t_Spike = vec(t_Spike)
  t_SpikeC = []
  t_SpikeC = vec(t_SpikeC)
  push!(t_Spike, freqE)
  push!(t_SpikeC, freqE+delay_pre)

  period = 0
  dayReset = 1


  # Step start and stop values
  # TstartE::Int64 = convert(Int64,TstepEinit/dt)
  # TstopE::Int64 = convert(Int64,TstepEfinal/dt)
  TstartI1::Int64 = convert(Int64,tstepIinit1/dt) 
  TstartI2::Int64 = convert(Int64,tstepIinit2/dt)
  TstartI3::Int64 = convert(Int64,tstepIinit3/dt)
  TstartI4::Int64 = convert(Int64,tstepIinit4/dt)
  TstopI::Int64   = convert(Int64,tstepIfinal/dt)

  #Nous tructural plasticity quand T2-->T3 et T4-->Fin

  for z = 1:Tdt

    for j = 1:ncells

      ## Excitatory cell
      if (j<=nEcells)
        Iapp = IappE
        if (z >= TstartI1 && z< TstartI2) || (z >= TstartI3 && z< TstartI4) ##EVEIL
          if mod(z*dt, freq) == 0
            newFreqE = rand(X)
            period += freq
            if (z >= TstartI1 && z< TstartI2)
              push!(t_Spike, newFreqE+period)
            end

            if  (z >= TstartI3 && z< TstartI4)
              if dayReset == 1
                period = 0
                dayReset = 0
              end
              push!(t_Spike, newFreqE+period+tstepIinit3)
            end
          end

          if length(t_Spike) > 0

            for i= 1:length(t_Spike)

              if round(t_Spike[i], digits = 2) == round(z*dt, digits = 2)|| spikedE == 1
                Iappstep = IstepE
                spikedE = 1
                if round(t_Spike[i] + duration, digits =2) == round(z*dt, digits = 2)
                  spikedE = 0
                  deleteat!(t_Spike, i)
                  Iappstep = 0.
                break
              end
            else
              Iappstep = 0.
            end
          end
        else
          Iappstep = 0.
        end
      else
        Iappstep = IstepE2
      end

      V[j] += dV(gNavec_E[j], gKdvec_E[j], glvec_E[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_E[j], gHvec_E[j], gKCavec_E[j])
      Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_E[j],k1vec_E[j], k2vec_E[j])
    end


      ## Inhibitory cell
    if (j > nEcells && j <= nEcells+nIcells)
      Iapp = IappI
      if z >= TstartI1 && z< TstartI2
        Iappstep = IstepI1

      elseif z >= TstartI2 && z< TstartI3
        Iappstep = IstepI2
        period = 0
      elseif z >= TstartI3 && z< TstartI4
        Iappstep = IstepI3
      elseif z >= TstartI4 && z< TstopI
        period = 0
        Iappstep = IstepI4
      else
        Iappstep = 0.
      end
        V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
    end

## Cortical cell
              if ( j > 1+nIcells && j <=ncells)
                Iapp = IappC

                if (z >= TstartI1 && z< TstartI2) || (z >= TstartI3 && z< TstartI4)
                  if mod(z*dt, freq) == 0
                    newFreqC = rand(X)
                    period += freq

                    if (z >= TstartI1 && z< TstartI2)
                      push!(t_SpikeC, newFreqC+period)
                    end

                    if  (z >= TstartI3 && z< TstartI4)
                      push!(t_SpikeC, newFreqC+period+tstepIinit3)
                    end
                  end

                  if length(t_SpikeC) > 0

                    for i= 1:length(t_SpikeC)

                      if round(t_SpikeC[i], digits = 2) == round(z*dt, digits = 2) || spikedC == 1
                        Iappstep = IstepC
                        spikedC = 1
                        if round(t_SpikeC[i] + duration, digits =2) == round(z*dt, digits = 2)
                          spikedC = 0
                          deleteat!(t_SpikeC, i)
                          Iappstep = 0.
                          break
                        end
                      else
                        Iappstep = 0.
                      end
                    end
                    else
                      Iappstep = 0.
                    end
                  else
                    Iappstep = IstepC2
                end


                V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
                Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
              end

## Interactions: excitation
if j == 1
    V[j] += (dt)*(1/C)*(-gIEGABAA*GABAA[2]*(Vprev[j]+70))
    V[j] += (dt)*(1/C)*(-gIEGABAB*GABAB[2]*(Vprev[j]+85))

end

if j == 3
    V[j] += (dt)*(1/C)*(-gEC*AMPA[1]*(Vprev[j]-0))
    V[j] += (dt)*(1/C)*(-gICGABAA*GABAA[2]*(Vprev[j]+70))
    V[j] += (dt)*(1/C)*(-gICGABAB*GABAB[2]*(Vprev[j]+85))
end

## Plasticity
      if (j <=nEcells)
        x_prev = x
        x += dx(x, V[j], Vprev[j])
        x2_prev = x2
        x2 += dx2(x2, V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          last = 1  #Pré en dernier
          count_post += 1
        end
      end

      if (j == 3)
        y_prev = y
        y += dy(y, V[j], Vprev[j])
        y2_prev = y2
        y2 += dy2(y2, V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          last = 0 #Post en dernier
          count_pre += 1
        end
      end
      if last == 1
        #Soft bound
        #new = A_w(w, Ap)*x+w
        # Hard bound
        new = -y_prev*w*(A2_m + A3_m*x2_prev) + w
        w_dot = -y_prev*w*(A2_m + A3_m*x2_prev) # for late-weight calculation
        # new = -y*(A2_m + A3_m*x2_prev) + w

        #new = -y*A2_m  + w

        w += dw(new, w, tau)

        #Calculation if the late-weight
        if (z >= TstartI2 && z< TstartI3)
          
          l += dl(w_dot, tau_l)
          
        
        elseif (z>= TstartI4) 
          l += dl(w_dot, tau_l)
        elseif (z >= TstartI1 && z< TstartI2) && (scenario == 2) #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        elseif (z >= TstartI3 && z< TstartI4) && (scenario == 2)  #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        else
          l += (dt)*0
        end

        



        # if (w> wMax)
        #   w= wMax
        # end
        #HB
        #A = A_w(w, Ap)

      elseif last == 0
        # Soft bound
        #new =  A_w(w, -Am)*y+w
        #Hard bound
        new = x_prev*(1-w)*(A2_p+ A3_p*y2_prev) + w
        w_dot = x_prev*(1-w)*(A2_p+ A3_p*y2_prev)

        # new = x*(A2_p+ A3_p*y2_prev) + w

        w += dw(new, w, tau)

        if (z >= TstartI2 && z< TstartI3)
          l += dl(w_dot, tau_l)
        elseif (z>= TstartI4) 
          l += dl(w_dot, tau_l)
        elseif (z >= TstartI1 && z< TstartI2) && (scenario == 2) #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        elseif (z >= TstartI3 && z< TstartI4) && (scenario == 2) #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        else
          l += (dt)*0
        end

        # if (w > wMax)
        #   w= wMax
        # end
       

      elseif last == 2
        w += dw(new, w, tau)

      
        
        #Hard Bound
        # if (w > wMax)
        #   w= wMax
        # end
      end
      last = 2


      gEC = gEC_const*w
      mNa[j] += dmNa(Vprev[j],mNa[j])
      hNa[j] += dhNa(Vprev[j],hNa[j])
      mKd[j] += dmKd(Vprev[j],mKd[j])
      mCaT[j] += dmCaT(Vprev[j],mCaT[j])
      hCaT[j] += dhCaT(Vprev[j],hCaT[j])
      mH[j] += dmH(Vprev[j],mH[j])

      AMPA[j] += dAMPA(Vprev[j],AMPA[j])
      GABAA[j] +=dGABAA(Vprev[j],GABAA[j])
      GABAB[j] +=dGABAB(Vprev[j],GABAB[j])

      Vprev = copy(V)
      
    end
    # VV[z,:]  = copy(V')
    if(mod(z,1/dt)==0)
      ww[idx_w] = w

      VV[idx_w,:]  = copy(V')
      ll[idx_w] = l

    

   



      idx_w+=1
    end

   
    

  end
 
  
  return ww, VV, ll
end




function simulate_scenario_conso(freq::Float64, delay_pre::Float64, w::Float64, l::Float64)


    gEC = gEC_const*w

    # Synaptic connections
    #SYN = [(gEE/4)*(rand(nEcells,nEcells)-0.5*ones(nEcells, nEcells)).+gEE (gEI/4)*(rand(nEcells,nIcells)-0.5*ones(nEcells, nIcells)).+gEI;(gIE/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE (gII/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII;(gIE2/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE2 (gII2/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII2]
    SYN = [gEE*ones(nEcells,nEcells) gEI*ones(nEcells,nIcells) gEC*ones(nEcells, nCcells);gIEGABAA*ones(nIcells,nEcells) gIIGABAA*ones(nIcells,nIcells) gICGABAA*ones(nIcells, nCcells);gIEGABAB*ones(nIcells,nEcells) gIIGABAB*ones(nIcells,nIcells) gICGABAB*ones(nIcells, nCcells);gCE*ones(nCcells,nEcells) gCI*ones(nCcells,nIcells) gCC*ones(nCcells,nCcells)]
    print(SYN[1,3],' ' , SYN[2,3])
    # Initial conditions
    V=-60*ones(ncells)
    Vprev=-60*ones(ncells)
    x::Float64 = 0.
    y::Float64 = 0.
    y2::Float64 = 0.
    x2::Float64 = 0.

    x_prev::Float64 = 0.
    y_prev::Float64 = 0.
    ww=zeros(convert(Int64,Tdt*dt))
    VV=zeros(convert(Int64,Tdt*dt), ncells)
    ll=0.1*ones(convert(Int64,Tdt*dt))

    idx_w=1
    # xx = zeros(Tdt)
    # yy = zeros(Tdt)
    # VV = zeros(Tdt,ncells)

    stop = 0

    mNa=mNainf(V[1])*ones(ncells)
    hNa=hNainf(V[1])*ones(ncells)
    mKd=mKdinf(V[1])*ones(ncells)
    mCaT=mCaTinf(V[1])*ones(ncells)
    hCaT=hCaTinf(V[1])*ones(ncells)
    mH=mHinf(V[1])*ones(ncells)
    Ca=((-k1_E/k2_E)*ICaT(V[1],mCaT[1],hCaT[1],gCaT_E[1]))*ones(ncells)
    AMPA=zeros(ncells)
    GABAA = zeros(ncells)
    GABAB = zeros(ncells)
    last = 2
    steady = 0.
    A = 0.
    new = w
    tau = 0.1
  

    w_dot=0

    zdt_freqE = 0
    z_delayed = delay_pre
    X = Distributions.Gaussian(delay_pre,0.1*delay_pre)
    delayC = rand(X)
    spikedE = 0
    spikedC = 0
    newFreqFlag = 0
    newfreqE = 0

    burst::Int64 = 0
    t_Spike = []
    t_Spike = vec(t_Spike)
    t_SpikeC = []
    t_SpikeC = vec(t_SpikeC)
    # push!(t_Spike, freq)
    # push!(t_SpikeC, freq+delayC)

    period = 0.
    dayReset = 1



    # Step start and stop values
    TstartI1::Int64 = convert(Int64,tstepIinit1/dt)
    TstartI2::Int64 = convert(Int64,tstepIinit2/dt)
    TstartI3::Int64 = convert(Int64,tstepIinit3/dt)
    TstartI4::Int64 = convert(Int64,tstepIinit4/dt)
    TstopI::Int64   = convert(Int64,tstepIfinal/dt)

    for z = 1:Tdt

      for j = 1:ncells

      ## Excitatory cell
      if (j<=nEcells)
        Iapp = IappE
        if (z >= TstartI1 && z< TstartI2) || (z >= TstartI3 && z< TstartI4)
        # z_delayed = z*dt + delay_pre
         # if (mod(z_delayed, freq) > freq - duration && z_delayed >=0 )
         #   Iappstep = IstepE
         # else
         #   Iappstep = 0.
          # end
          if mod(z*dt, freq) == 0

             newDelay = rand(X)
             period += freq

            if (z >= TstartI1 && z< TstartI2)
              push!(t_Spike, period)
              push!(t_SpikeC, newDelay+period)
            end

        if  (z >= TstartI3 && z< TstartI4)
          if dayReset == 1
            period = 0
            dayReset =0
          end
          push!(t_Spike, period+tstepIinit3)
          push!(t_SpikeC, newDelay+period+tstepIinit3)
        end
      end

      if length(t_Spike) > 0
        for i= 1:length(t_Spike)
          if round(t_Spike[i], digits = 2) == z*dt || spikedE == 1

            Iappstep = IstepE
            spikedE = 1
            if round(t_Spike[i]+duration, digits = 2) == z*dt
              spikedE = 0
              deleteat!(t_Spike, i)
              Iappstep = 0.
              break
            end
          else
            Iappstep = 0.
          end
        end
      else
        Iappstep = 0.
      end


    else
      Iappstep = IstepE2
    end

    V[j] += dV(gNavec_E[j], gKdvec_E[j], glvec_E[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_E[j], gHvec_E[j], gKCavec_E[j])
    Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_E[j],k1vec_E[j], k2vec_E[j])
  end


  ## Inhibitory cell
  if (j > nEcells && j <= nEcells+nIcells)
    Iapp = IappI
    if z >= TstartI1 && z< TstartI2
      Iappstep = IstepI1
    elseif z >= TstartI2 && z< TstartI3
      Iappstep = IstepI2
    elseif z >= TstartI3 && z< TstartI4
      Iappstep = IstepI3
    elseif z >= TstartI4 && z< TstopI
      Iappstep = IstepI4
    else
      Iappstep = 0.
    end
    V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
    Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
  end

  ## Cortical cell
  if ( j > 1+nIcells && j <=ncells)
    Iapp = IappC
    if (z >= TstartI1 && z< TstartI2) || (z >= TstartI3 && z< TstartI4)

      # z_delayed = z*dt
      # if (mod(z_delayed, freq) > freq - duration && z_delayed >=0 )
      #   Iappstep = IstepC
      # else
      #   Iappstep = 0.
      # end
      if length(t_SpikeC) >0
        for i = 1:length(t_SpikeC)
          if round(t_SpikeC[i], digits = 2) == round(z*dt, digits = 2) || spikedC == 1
            Iappstep = IstepC
            spikedC = 1
            if round(t_SpikeC[i] + duration, digits = 2) == round(z*dt, digits = 2)
              spikedC = 0
              deleteat!(t_SpikeC, i)
              Iappstep = 0.
              break
            end
          else
            Iappstep = 0.
          end
        end
      else
        Iappstep = 0.

      end

    else
      Iappstep = IstepC2
    end

    V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
    Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
    end


      ## Interactions: excitation
    if j == 1
      V[j] += (dt)*(1/C)*(-gIEGABAA*GABAA[2]*(Vprev[j]+70))
      V[j] += (dt)*(1/C)*(-gIEGABAB*GABAB[2]*(Vprev[j]+85))

   end

    if j == 3
      V[j] += (dt)*(1/C)*(-gEC*AMPA[1]*(Vprev[j]-0))
      V[j] += (dt)*(1/C)*(-gICGABAA*GABAA[2]*(Vprev[j]+70))
      V[j] += (dt)*(1/C)*(-gICGABAB*GABAB[2]*(Vprev[j]+85))
    end


        ## Plasticity
      if (j <=nEcells)
        x_prev = x
        x += dx(x, V[j], Vprev[j])
        x2_prev = x2
        x2 += dx2(x2, V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          last = 1  #Pre in last
        end
      end

      if (j == 3)
        y_prev = y
        y += dy(y, V[j], Vprev[j])
        y2_prev = y2
        y2 += dy2(y2, V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          last = 0 #Post in last
        end
      end

      if last == 1 #Pre in last
        #Soft bound
        #new = A_w(w, Ap)*x+w
        # Hard bound
        new = -y_prev*w*(A2_m + A3_m*x2_prev) + w
        w_dot = -y_prev*w*(A2_m + A3_m*x2_prev)
        # new = -y*(A2_m + A3_m*x2_prev) + w

        #new = -y*A2_m  + w

        w += dw(new, w, tau)

        # Calculation of the late-weight
        if (z >= TstartI2 && z< TstartI3)
          l += dl(w_dot, tau_l)
        elseif (z>= TstartI4) 
          l += dl(w_dot, tau_l)
        elseif (z >= TstartI1 && z< TstartI2) && (scenario == 2) #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        elseif (z >= TstartI3 && z< TstartI4) && (scenario == 2) #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        else
          l += (dt)*0
        end
       
      

        # if (w> wMax)
        #   w= wMax
        # end
        #HB
        #A = A_w(w, Ap)

      elseif last == 0 #Post in last
        # Soft bound
        #new =  A_w(w, -Am)*y+w
        #Hard bound
        new = x_prev*(1-w)*(A2_p+ A3_p*y2_prev) + w
        w_dot = x_prev*(1-w)*(A2_p+ A3_p*y2_prev)

        # new = x*(A2_p+ A3_p*y2_prev) + w

        w += dw(new, w, tau)

        if (z >= TstartI2 && z< TstartI3)
         
        
          l += dl(w_dot, tau_l)
        
        elseif (z>= TstartI4) 
          l += dl(w_dot, tau_l)
        elseif (z >= TstartI1 && z< TstartI2) && (scenario == 2) #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        elseif (z >= TstartI3 && z< TstartI4) && (scenario == 2) #when structural plasticity is also implemented in tonic
          l += dl_tonic(w_dot, tau_l)
        else
          l += (dt)*0
        end

        

        # if (w > wMax)
        #   w= wMax
        # end

      elseif last == 2
        w += dw(new, w, tau)

        
        #Hard Bound
        # if (w > wMax)
        #   w= wMax
        # end
      end
      last = 2

      gEC = gEC_const*w
      mNa[j] += dmNa(Vprev[j],mNa[j])
      hNa[j] += dhNa(Vprev[j],hNa[j])
      mKd[j] += dmKd(Vprev[j],mKd[j])
      mCaT[j] += dmCaT(Vprev[j],mCaT[j])
      hCaT[j] += dhCaT(Vprev[j],hCaT[j])
      mH[j] += dmH(Vprev[j],mH[j])

      AMPA[j] += dAMPA(Vprev[j],AMPA[j])
      GABAA[j] +=dGABAA(Vprev[j],GABAA[j])
      GABAB[j] +=dGABAB(Vprev[j],GABAB[j])

      Vprev = copy(V)
      
    end
    # VV[z,:]  = copy(V')
    



    if(mod(z,1/dt)==0)
      ww[idx_w] = w
      VV[idx_w, :] = copy(V')
      ll[idx_w] = l

    


      

      idx_w+=1

    end
    


  
    # ww[z] = copy(w)
  end

  #println("Max w = ", maximum(ww))
  #println("Max g = ", maximum(ll))

  return ww, VV, ll
end


#--#

function debug_scenario_conso(freq::Float64, delay_pre::Float64, w::Float64)


    gEC = gEC_const*w

    # Synaptic connections
    #SYN = [(gEE/4)*(rand(nEcells,nEcells)-0.5*ones(nEcells, nEcells)).+gEE (gEI/4)*(rand(nEcells,nIcells)-0.5*ones(nEcells, nIcells)).+gEI;(gIE/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE (gII/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII;(gIE2/4)*(rand(nIcells,nEcells)-0.5*ones(nIcells, nEcells)).+gIE2 (gII2/4)*(rand(nIcells,nIcells)-0.5*ones(nIcells, nIcells)).+gII2]
    SYN = [gEE*ones(nEcells,nEcells) gEI*ones(nEcells,nIcells) gEC*ones(nEcells, nCcells);gIEGABAA*ones(nIcells,nEcells) gIIGABAA*ones(nIcells,nIcells) gICGABAA*ones(nIcells, nCcells);gIEGABAB*ones(nIcells,nEcells) gIIGABAB*ones(nIcells,nIcells) gICGABAB*ones(nIcells, nCcells);gCE*ones(nCcells,nEcells) gCI*ones(nCcells,nIcells) gCC*ones(nCcells,nCcells)]
    print(SYN[1,3],' ' , SYN[2,3])
    # Initial conditions
    V=-60*ones(ncells)
    Vprev=-60*ones(ncells)
    x::Float64 = 0.
    y::Float64 = 0.
    y2::Float64 = 0.
    x2::Float64 = 0.

    x_prev::Float64 = 0.
    y_prev::Float64 = 0.
    ww=zeros(Tdt)
    idx_w=1
    # xx = zeros(Tdt)
    # yy = zeros(Tdt)

    stop = 0

    mNa=mNainf(V[1])*ones(ncells)
    hNa=hNainf(V[1])*ones(ncells)
    mKd=mKdinf(V[1])*ones(ncells)
    mCaT=mCaTinf(V[1])*ones(ncells)
    hCaT=hCaTinf(V[1])*ones(ncells)
    mH=mHinf(V[1])*ones(ncells)
    Ca=((-k1_E/k2_E)*ICaT(V[1],mCaT[1],hCaT[1],gCaT_E[1]))*ones(ncells)
    AMPA=zeros(ncells)
    GABAA = zeros(ncells)
    GABAB = zeros(ncells)
    last = 2
    steady = 0.
    A = 0.
    new = w
    tau = 0.1


    zdt_freqE = 0
    z_delayed = delay_pre
    X = Distributions.Gaussian(freq,0.3*freq)
    freqE = rand(X)
    spikedE = 0
    spikedC = 0
    newFreqFlag = 0
    newfreqE = 0

    burst::Int64 = 0
    t_Spike = []
    t_Spike = vec(t_Spike)
    t_SpikeC = []
    t_SpikeC = vec(t_SpikeC)
    # push!(t_Spike, freq)
    # push!(t_SpikeC, freq+delayC)

    period = 0.
    dayReset = 1


    # Step start and stop values
    TstartI1::Int64 = convert(Int64,tstepIinit1/dt)
    TstartI2::Int64 = convert(Int64,tstepIinit2/dt)
    TstartI3::Int64 = convert(Int64,tstepIinit3/dt)
    TstartI4::Int64 = convert(Int64,tstepIinit4/dt)
    TstopI::Int64   = convert(Int64,tstepIfinal/dt)
    for z = 1:Tdt

      for j = 1:ncells

  ## Excitatory cell
        if (j<=nEcells)
          Iapp = IappE
            if (z >= TstartI1 && z< TstartI2) || (z >= TstartI3 && z< TstartI4)
              if mod(z*dt, freq) == 0
                newFreqE = rand(X)
                period += freq
                if (z >= TstartI1 && z< TstartI2)
                  push!(t_Spike, newFreqE+period)
                end

                if  (z >= TstartI3 && z< TstartI4)
                  if dayReset == 1
                    period = 0
                    dayReset = 0
                  end
                  push!(t_Spike, newFreqE+period+tstepIinit3)
                end
              end

              if length(t_Spike) > 0

                for i= 1:length(t_Spike)

                  if round(t_Spike[i], digits = 2) == round(z*dt, digits = 2)|| spikedE == 1
                    Iappstep = IstepE
                    spikedE = 1
                    if round(t_Spike[i] + duration, digits =2) == round(z*dt, digits = 2)
                      spikedE = 0
                      deleteat!(t_Spike, i)
                      Iappstep = 0.
                      break
                    end
                  else
                    Iappstep = 0.
                  end
                end
                else
                  Iappstep = 0.
                end
              else
                Iappstep = IstepE2
            end

            V[j] += dV(gNavec_E[j], gKdvec_E[j], glvec_E[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_E[j], gHvec_E[j], gKCavec_E[j])
            Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_E[j],k1vec_E[j], k2vec_E[j])
          end


  ## Inhibitory cell
        if (j > nEcells && j <= nEcells+nIcells)
          Iapp = IappI
          if z >= TstartI1 && z< TstartI2
            Iappstep = IstepI1

          elseif z >= TstartI2 && z< TstartI3
            Iappstep = IstepI2
            period = 0
          elseif z >= TstartI3 && z< TstartI4
            Iappstep = IstepI3
          elseif z >= TstartI4 && z< TstopI
            period = 0
            Iappstep = IstepI4
          else
            Iappstep = 0.
          end
          V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
          Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
        end

  ## Cortical cell
        if ( j > 1+nIcells && j <=ncells)
          Iapp = IappC

          if (z >= TstartI1 && z< TstartI2) || (z >= TstartI3 && z< TstartI4)
            if mod(z*dt, freq) == 0
              newFreqC = rand(X)
              period += freq

              if (z >= TstartI1 && z< TstartI2)
                push!(t_SpikeC, newFreqC+period)
              end

              if  (z >= TstartI3 && z< TstartI4)
                push!(t_SpikeC, newFreqC+period+tstepIinit3)
              end
            end

            if length(t_SpikeC) > 0

              for i= 1:length(t_SpikeC)

                if round(t_SpikeC[i], digits = 2) == round(z*dt, digits = 2) || spikedC == 1
                  Iappstep = IstepC
                  spikedC = 1
                  if round(t_SpikeC[i] + duration, digits =2) == round(z*dt, digits = 2)
                    spikedC = 0
                    deleteat!(t_SpikeC, i)
                    Iappstep = 0.
                    break
                  end
                else
                  Iappstep = 0.
                end
              end
              else
                Iappstep = 0.
              end
            else
              Iappstep = IstepC2
          end


          V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
          Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
        end



  ## Interactions: excitation
  if j == 1
      V[j] += (dt)*(1/C)*(-gIEGABAA*GABAA[2]*(Vprev[j]+70))
      V[j] += (dt)*(1/C)*(-gIEGABAB*GABAB[2]*(Vprev[j]+85))

  end

  if j == 3
      V[j] += (dt)*(1/C)*(-gEC*AMPA[1]*(Vprev[j]-0))
      V[j] += (dt)*(1/C)*(-gICGABAA*GABAA[2]*(Vprev[j]+70))
      V[j] += (dt)*(1/C)*(-gICGABAB*GABAB[2]*(Vprev[j]+85))
  end


## Plasticity
      if (j <=nEcells)
        x_prev = x
        x += dx(x, V[j], Vprev[j])
        x2_prev = x2
        x2 += dx2(x2, V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          last = 1  #Pre in last
        end
      end

      if (j == 3)
        y_prev = y
        y += dy(y, V[j], Vprev[j])
        y2_prev = y2
        y2 += dy2(y2, V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          last = 0 #Post in last
        end
      end
      if last == 1
        #Soft bound
        #new = A_w(w, Ap)*x+w
        # Hard bound
        new = -y_prev*w*(A2_m + A3_m*x2_prev) + w
        # new = -y*(A2_m + A3_m*x2_prev) + w

        #new = -y*A2_m  + w

        w += dw(new, w, tau)

        # if (w> wMax)
        #   w= wMax
        # end
        #HB
        #A = A_w(w, Ap)

      elseif last == 0
        # Soft bound
        #new =  A_w(w, -Am)*y+w
        #Hard bound
        new = x_prev*(1-w)*(A2_p+ A3_p*y2_prev) + w

        # new = x*(A2_p+ A3_p*y2_prev) + w

        w += dw(new, w, tau)

        # if (w > wMax)
        #   w= wMax
        # end

      elseif last == 2
        w += dw(new, w, tau)
        #Hard Bound
        # if (w > wMax)
        #   w= wMax
        # end
      end
      last = 2

      gEC = gEC_const*w
      mNa[j] += dmNa(Vprev[j],mNa[j])
      hNa[j] += dhNa(Vprev[j],hNa[j])
      mKd[j] += dmKd(Vprev[j],mKd[j])
      mCaT[j] += dmCaT(Vprev[j],mCaT[j])
      hCaT[j] += dhCaT(Vprev[j],hCaT[j])
      mH[j] += dmH(Vprev[j],mH[j])

      AMPA[j] += dAMPA(Vprev[j],AMPA[j])
      GABAA[j] +=dGABAA(Vprev[j],GABAA[j])
      GABAB[j] +=dGABAB(Vprev[j],GABAB[j])

      Vprev = copy(V)
    end
    VV[z,:]  = copy(V')
    # if(mod(z,1/dt)==0)
    #   ww[idx_w] = w
    #   idx_w+=1
    # end
    ww[z] = copy(w')
  end
  return VV, ww
end
