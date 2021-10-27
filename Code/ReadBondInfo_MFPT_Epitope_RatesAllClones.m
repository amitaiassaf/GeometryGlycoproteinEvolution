%  Author:  Assaf Amitai
%            amitaiassaf@gmail.com
%            https://github.com/amitaiassaf/SpikeGeometry


% The script computes the on-rates for the first and second antibody (Ab)
% arm binding to the respective epitope using output (dump files) from
% Lammps runs. The script further computes the probability of such
% events.

% To run the script put all the dump files in the run directory of the
% script. The dump file should have the form:
%  SpikeSim_idum_${idum}_virushaSpikes${SpikeN}Bond_Ep${epitopenum}_Cond${Cond}
% ${idum}: seed number for the random number generator of Lammps.
% ${SpikeN}: the number of spikes (HA molecules or S proteins) on the surface of the virus.
% ${epitope}: the target residue (epitope) on the surface of HA. There is a total of 228 epitopes for HA
% ${Cond}: a variable indicating different initial positions of the Ab.


% The output of the script is written to a matflie which contains the following variables
% P1Arm: the probability of one Ab arm binding to the epitope
% P2Arm: the probability of a second Ab arm binding to the epitope, given that the first is already bound.
% kOn1Arm: the on-rate for the first arm binding.
% kOn2Arm: the on-rate for the second arm binding, given that the first is already bound.


dt = 1e-7; % conversion to the units of time. Based on Delta t chosen in the MD simulations

SpikeN = 40;

CondNum = 7; % The number of different initial configurations of the Ab with respect to the immunogen

idumArr = [400:402] % idum span (seed to the random number generator of the simulation)

epitopeArr = [7:234]; % All possible epitopes on HA (total of 228 epitopes)

MFET_1stArm_Arr = zeros(1,length(epitopeArr)); % the mean first encounter time of 1 arm of the Ab to bind to its epitope
MFET_2ndArm_Arr = zeros(1,length(epitopeArr)); % the mean first encounter time of both arms of the Ab to bind to their epitope
FET_1stArm_ArrCell = cell(1,length(epitopeArr)); % cell array containing all the first encounter times of 1 arm of the Ab to bind to its epitope
FET_2ndArm_ArrCell = cell(1,length(epitopeArr)); % cell array containing all the first encounter times of both arsm of the Ab to bind to their epitope
SimNumEp = zeros(1,length(epitopeArr)); % number of simulations for each epitope

for ep=1:length(epitopeArr)
    
    epitope = epitopeArr(ep);

    FET_1stArm_epitope = [];
    FET_2ndArm_epitope = [];

    for cond=1:CondNum
        
        for idum =idumArr
            TEF_1stArm_i = [];
            TEF_2ndArm_i = [];
            
            filename = ['SpikeSim_idum_',num2str(idum),'_virushaSpikes',num2str(SpikeN),'_Bond_Ep_',num2str(epitope),'_Cond',num2str(cond)];
            filename
            
            if(exist(filename)==0)
                continue;
            end
            
            SimNumEp(1,ep) = SimNumEp(1,ep) + 1;
            LatestArrvialtime = 0;
            Prev1ArmMFPT = 0;
            
            data = readdump_all_bonds(filename);
            
            OneArmBindTime = [];
            if(~isempty(find(data.BondNum)))
                Bondsidx = find(data.BondNum);
                
                Target1 = [5 epitope];
                Target2 = [epitope 5];
                
                for i=1:length(Bondsidx)
                    bonds = data.Bonds{Bondsidx(i)};
                    idx5 = find( (bonds(:,4)==5) | (bonds(:,5)==5) ); %% All attachments of Ab to epitope
                    if(length(idx5) & isempty(TEF_1stArm_i))
                        FET = dt*data.timestep(Bondsidx(i)); % The first encounter time of 1 Ab arm to target epitope
                        TEF_1stArm_i =  FET;
                        TEF_2ndArm_i = NaN;
                    end
                    
                    Arm1idx = 0;
                    Arm2idx = 0;
                    SpikeArm1 = 0;
                    SpikeArm2 = 0;
                    
                    OneArmBindTime = [OneArmBindTime ; data.timestep(Bondsidx(i)) ];
                    for b=1:length(idx5)
                        bond = bonds(idx5(b),:);
                        
                        if( isequal(bond(4:5),Target1) | isequal(bond(4:5),Target2) );

                            Position = find(bond(4:5)==5);
                            EpitopePosition = find(bond(4:5)==epitope);
                            idxArm = bond(Position+1);
                            idxSpike = bond(EpitopePosition+1);
                            if(Arm1idx==0)
                                Arm1idx = idxArm;
                                SpikeArm1 = idxSpike;
                            elseif(idxArm~=Arm1idx)
                                
                                if(isnan(TEF_2ndArm_i))
                                    
                            
                                    
                                    d1ArmBound = diff(OneArmBindTime);
                                    DTSim = min(d1ArmBound);
                                    dtTime = 1.05*min(d1ArmBound);
                                    d1ArmArrival = find(d1ArmBound>dtTime);
                                    if(length(d1ArmArrival))
                                        LatestArrvialtime = OneArmBindTime(d1ArmArrival(end)+1);
                                    else
                                        LatestArrvialtime = FET;
                                    end
                                    
                                    FET2arm = dt*(data.timestep(Bondsidx(i))-LatestArrvialtime); % the time it takes the second Ab arm to bind once the first binds
                                    if(FET2arm==0) % If the second arm binds at the same simlation step as the first arm impute the 2nd arm binding bind time as dt
                                        FET2arm = dt*DTSim;
                                        TEF_2ndArm_i = FET2arm;
                                    else
                                        if(Prev1ArmMFPT~=LatestArrvialtime)
                                            if(isnan(TEF_2ndArm_i))
                                                TEF_2ndArm_i = FET2arm;
                                            else
                                                TEF_2ndArm_i = [TEF_2ndArm_i ; FET2arm];
                                            end
                                        end
                                        
                                    end
                                    Prev1ArmMFPT = LatestArrvialtime;
                                end

                            end
         
                        end
                    end
                end
            end
            FET_1stArm_epitope = [FET_1stArm_epitope ; TEF_1stArm_i];
            FET_2ndArm_epitope = [FET_2ndArm_epitope ; TEF_2ndArm_i];
            
        end
    end
    
    MFET_1stArm_Arr(1,ep) = mean(FET_1stArm_epitope);
    idxnotnan = find(~isnan(FET_2ndArm_epitope));
    MFET_2ndArm_Arr(1,ep) = mean(FET_2ndArm_epitope(idxnotnan));
    
    FET_1stArm_ArrCell{1,ep} = FET_1stArm_epitope;
    FET_2ndArm_ArrCell{1,ep} = FET_2ndArm_epitope;
    
end


P1Arm = zeros(1,length(epitopeArr)); % The fraction of simulations where 1 arm binds
P2Arm = zeros(1,length(epitopeArr));% The fraction of simulations where 2 arms bind

for j=1:length(epitopeArr)
    OneArmMFPT = FET_1stArm_ArrCell{1,j};
    idxnotnan = find(~isnan(OneArmMFPT));
    P1Arm(1,j) = length(idxnotnan)/SimNumEp(1,j);
    
    TwoArmMFPT = FET_2ndArm_ArrCell{1,j};
    idxnotnan = find(~isnan(TwoArmMFPT));
    P2Arm(1,j) = length(idxnotnan)/SimNumEp(1,j);    
end

kOn1Arm = zeros(1,length(epitopeArr)); % The on-rate of the first Ab arm
kOn2Arm = zeros(1,length(epitopeArr)); % The on-rate of the second arm of the Ab once the first arm binds

for j=1:length(epitopeArr)
    OneArmMFPT = FET_1stArm_ArrCell{1,j};
    idxnotnan = find(~isnan(OneArmMFPT));
    if(length(idxnotnan))
        kOn1Arm(1,j) = 1/mean(OneArmMFPT(idxnotnan));
    end
    
    TwoArmMFPT = FET_2ndArm_ArrCell{1,j};
    idxnotnan = find(~isnan(TwoArmMFPT));
    if(length(idxnotnan))
        kOn2Arm(1,j) = 1/mean(TwoArmMFPT(idxnotnan));
    end
end
 
flag_str = ['virusha',num2str(SpikeN),'Spikes'];

save([flag_str,'.mat'],'P1Arm','P2Arm','kOn1Arm','kOn2Arm','-v7.3');
