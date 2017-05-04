%%MC-1 v3 - BY ANDY OGILVY & PAUL IONELE

%Consider only photon transport. All charged particles produced deposit
%their kinetic energy at their point of origin.

%Refer to report Problem 1 section for simulation information.

clear; clc

%Called function returns arrays of energies and atten. coeffs for various
%interactions. Note a single unique index specifies cooresponding values.
file1 = 'mass_atten.txt';
[energies,coherent,comptons,photoels,pairtrip,energytr,energyab,radifrac] = table1(file1);

%Expressing mass atten. coeffs. as relative probabilites (for plotting).
coherent_p = zeros(1,length(energies));
comptons_p = zeros(1,length(energies));
photoels_p = zeros(1,length(energies));
pairtrip_p = zeros(1,length(energies));
for k = 1:length(energies)
    tma1 = coherent(k) + comptons(k) + photoels(k) + pairtrip(k);
    coherent_p(k) = coherent(k)/tma1;
    comptons_p(k) = comptons(k)/tma1;
    photoels_p(k) = photoels(k)/tma1;
    pairtrip_p(k) = pairtrip(k)/tma1;
end

%Arbitrary
N = 10; %number of particles
ec = 0.0001; %energy cutoff in MeV (can/should set to 0.01)
theta = 0:0.01:pi; %range of angles between [0,180]
times = zeros(1,length(energies)); %to record computation time

%Counters, as an array because we will have different #s for each energy.
%A given photon may interact more than one (i.e. compton or elastic).
%Preallocate for speed and because reasons.
ens = energies; %new variable. Was supposed to be in (10keV) increments.
coher_count = zeros(1,length(ens));
compt_count = zeros(1,length(ens));
photo_count = zeros(1,length(ens));
pairp_count = zeros(1,length(ens));

compt_e_trans = zeros(1,length(ens)); %energy transferred to electrons
num_interaction = zeros(1,length(ens)); %total # ints at each energy
ticker = zeros(1,length(energies)); %for the avg. # of ints a photon undergoes

%The following array dramatically slows the simulation. Remove if not req.
compt_ang_dist = []; %angular distribution compton scattered photons

for i = 1:length(energies) %loop over all energies in list
    tic;
    e = energies(i);
    if e <= ec
        %No need to sub-loop over energies we don't want.
        continue
    end
    
    for j = 1:N %N photons to be simulated
        e = energies(i); %set new photon to a list energy
        while e >= ec %while jth photon energy is greater than cutoff energy 
            ticker(i) = ticker(i) + 1;

            %Linear interpolation between mass attenuation coefficients.
            coherent_interp = interp1(energies, coherent, e);
            comptons_interp = interp1(energies, comptons, e);
            photoels_interp = interp1(energies, photoels, e);
            pairtrip_interp = interp1(energies, pairtrip, e);
            
            %Total mass attenuation for e.
            tma = coherent_interp + comptons_interp +...
                photoels_interp + pairtrip_interp; 
            
            %Probabilities for interactions, used to select interaction.
            p1 = coherent_interp/tma;
            p2 = p1 + comptons_interp/tma;
            p3 = p2 + photoels_interp/tma;
            p4 = p3 + pairtrip_interp/tma;
           
            r1 = rand(1); %generate uniform random deviate on [0,1].
            
            if (r1 < p1)
                %COHERENT SCATTER; e_in = e_out, so nothing happens here!
                %Update counter for coherent scatter.
                
                [~,t_index] = min(abs(ens-e)); %index of closest energy in ens array
                coher_count(t_index) = coher_count(t_index) + 1; %+1 count increment at that energy
                num_interaction(t_index) = num_interaction(t_index) + 1;
                
            elseif (r1 >= p1) && (r1 < p2)
                %COMPTON SCATTER; scatter angle acc. to KN. This is the
                %only interaction in this simulation where the photon loses
                %energy and continues on.
                
                %Update counters for incident photon energy.
                [~,t_index] = min(abs(ens-e)); %index of closest energy
                compt_count(t_index) = compt_count(t_index) + 1; %count at that energy
                num_interaction(t_index) = num_interaction(t_index) + 1;
                
                %The following is to determine new photon energy. First, an
                %array of photon energies ep, based on photon scatter
                %angles [0,180] with fine sampling; theta defined in
                %header.
                ep = e./(1+(e/0.51099)*(1-cos(theta))); 
                
                %Here, we build the differential compton cross section
                %for the e MeV photon energy. This is again an array. We
                %will sample from this array.
                d_cr = (2*pi*sin(theta)).*((2.8179*10^-15)^2)/2.*((ep/e).^2).*(e./ep + ep/e - (sin(theta).^2));
                d_cr_max = max(d_cr); %returns maximum value of d_cr; no index
                %[~,c] = ismember(d_cr_max,d_cr); %call c for index of max
                
                while 1
                    %Using the 'rejection method' to determine energy.
                    %This loop will continue as long as r2 is rejected.
                    %Loop terminates by 'break'.
                    r2 = rand(1)*pi; %uniform random deviate on [0,pi]
                    r3 = rand(1)*d_cr_max; %[0,max_cross_section]

                    %Get scattered photon energy ep2 at the random angle
                    %r2. Evaluate d_cr again but for the random angle and
                    %cooresponding photon energy ep2 -> d_cr2.
                    ep2 = e/(1+(e/0.51099)*(1-cos(r2))); 
                    d_cr2 = (2*pi*sin(r2))*((2.8179*10^-15)^2)/2*((ep2/e)^2)*(e/ep2 + ep2/e - (sin(r2)^2));

                    if r3 <= (d_cr2)
                        %Accept r2 random angle, continue on to calculate
                        %the kinematic parameters of the scattered photon
                        %and Compton electron. Determine angles.
                                     
                        %Mean energy transferred to electron as a function
                        %of incident photon energy. Add to total fraction
                        %at an energy. Will then divide by total number of
                        %compton interactions at that energy.
                        t = e - ep2; %electron kinetic energy = initial photon e - scattered e
                        tf = t/e; %fraction of energy transferred to electron
                        [~,t_index] = min(abs(ens-e));
                        compt_e_trans(t_index) = compt_e_trans(t_index) + tf;

                        compt_ang_dist = [compt_ang_dist, r2];
                        
                        e = ep2; %update energy to new scattered photon energy

                        break %breaks out of while loop
                    else
                        %Reject r2 random angle, run while loop again.
                    end
                end           
            elseif (r1 >= p2) && (r1 < p3)
                %PHOTOELECTRIC; photon absorbed, electron ejected takes
                %energy according to e_elec = e - 543.1 eV.
                
                [~,t_index] = min(abs(ens-e));
                photo_count(t_index) = photo_count(t_index) + 1;
                num_interaction(t_index) = num_interaction(t_index) + 1;
                
                e = 0;                
            else
                %PAIR + TRIPLET PRODUCTION; photon disappears.
                
                [~,t_index] = min(abs(ens-e));
                pairp_count(t_index) = pairp_count(t_index) + 1;
                num_interaction(t_index) = num_interaction(t_index) + 1;
                
                e = 0;
            end
            %At this point, photon may have new energy dependent on the
            %interaction that it underwent.
        end
            %Deposit energy at current position. Nothing done here.
    end
    
    disp(energies(i)) %display progress
    times(i) = toc;
end

%Relative probabilities for interaction from simulation.
coher_rp = coher_count./num_interaction;
compt_rp = compt_count./num_interaction;
photo_rp = photo_count./num_interaction;
pairp_rp = pairp_count./num_interaction;

%Figures:
fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;
fig5 = figure;

%Generated plots.
figure(fig1);
loglog(ens,coher_rp,'-r','linewidth',1.5,'DisplayName', 'Coherent')
hold on
loglog(ens,compt_rp,':g','linewidth',1.5,'DisplayName', 'Compton')
loglog(ens,photo_rp,'-.b','linewidth',1.5,'DisplayName', 'Photoelectric')
loglog(ens,pairp_rp,'-k','linewidth',1.5,'DisplayName', 'Pair + Triple')
xlim([0,50])

%Actual plots.
loglog(energies,coherent_p,'DisplayName', 'Coherent - A')
loglog(energies,comptons_p,'DisplayName', 'Compton - A')
loglog(energies,photoels_p,'DisplayName', 'Photoelectric - A')
loglog(energies,pairtrip_p,'DisplayName', 'Pair + Triple - A')
legend('show','Location','South')
hold off

%Average fraction of energy transferred to recoil electrons in Compton
%interactions, as a function of the incident photon energy.
compt_e_trans_frac = compt_e_trans./compt_count;

figure(fig2);
semilogx(ens, compt_e_trans_frac)
xlabel('Energy [MeV]')
ylabel('Average Fraction of Energy Transferred')
title('Avg. Fraction of Energy Trans. to Electrons vs. Energy [Mev]')
axis([0,50,0,1])
grid('on')

%Average number of interactions a photon undergoes until it is absorbed in
%an infinite medium of water, as a funtion of photon energy.
figure(fig3);
semilogx(energies, ticker/N)
xlabel('Energy [MeV]')
ylabel('Average Number Interactions per Primary Incident Photon')
title('Average Number of Photon Interactions vs. Energy [Mev]')
axis([0,50,0,inf])
grid('on')

% Angular distribution of Compton scattered electrons.
figure(fig4);
hist(compt_ang_dist,100)
xlabel('Angle [radians]')
ylabel('Counts per Binned Angles')
title('Histogram for Angular Distribution of Compton Photons')
axis([0,pi,0, inf])

%Time per loop (if not interested comment out).
figure(fig5);
plot(times,'-ro')
xlabel('Energy/Loop Number')
ylabel('Time [s]')
title('Time [s] per Loop')
axis([0,39,0, inf])
grid('on')