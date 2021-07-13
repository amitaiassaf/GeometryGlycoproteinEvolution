%  Author :  Assaf Amitai
%            amitaiassaf@gmail.com
%            https://github.com/amitaiassaf/SpikeGeometry


dt = 1e-7;
BondsArr = [];
Boundscount = 0;
count = 1;
epitope = 11;

SpikesArr = [20];

NArr = [40];

epitopeArr = [7:234];

Boundscount = [];
BoundscountArr = zeros(length(NArr),length(epitopeArr));
OneArmcountArr = zeros(length(NArr),length(epitopeArr));
TwoArmcountArr = zeros(length(NArr),length(epitopeArr));

KonArr = zeros(length(NArr),length(epitopeArr));
Kon2ndArmArr = zeros(length(NArr),length(epitopeArr));
KonArrCell = cell(length(NArr),length(epitopeArr));
Kon2ndArmArrCell = cell(length(NArr),length(epitopeArr));

BoundscountArrNorm = zeros(length(NArr),length(epitopeArr));
BoundsTotcountArr = zeros(length(NArr),length(epitopeArr));
BondEnergytArr = zeros(length(NArr),length(epitopeArr));
BondDistArr = zeros(length(NArr),length(epitopeArr));

TwoArmDistBoundArr = zeros(length(NArr),length(epitopeArr));
TwoArmDistBoundSameSpikeArr = zeros(length(NArr),length(epitopeArr));
for n=1:length(NArr)
    N = NArr(n);

    SpikeNum = NArr(n);
    for ep=1:length(epitopeArr)
        
        epitope = epitopeArr(ep);
        Boundscount = 0;
        OneArmcount = 0;
        TwoArmcount = 0;
        TwoArmDistBoundSameSpike = 0;
        KonN_epitope = [];
        KonN_2ndArm_epitope = [];
        BondEnergy = [];
        BondDist = [];
        TwoArmDistBound = [];
        for cond=1:7

            for idum =[400:402]
                Kon = [];
                Kon2ndarm = [];
                

                filename = ['SpikeSim_idum_',num2str(idum),'_virushaSpikes',num2str(SpikeNum),'_Bond_Ep_',num2str(epitope),'_Cond',num2str(cond)];
                filename

                if(exist(filename)==0)
                    continue;
                end
                BoundscountArrNorm(n,ep) = BoundscountArrNorm(n,ep) + 1;
                LatestArrvialtime = 0;
                Prev1ArmMFPT = 0;
                
                data = readdump_all_bonds(filename);
                BoundsTotcountArr(n,ep) = BoundsTotcountArr(n,ep) + length(data.timestep);
            
                LatestMFPT = 0;
                OneArmBindTime = [];
                if(~isempty(find(data.BondNum)))
                    Bondsidx = find(data.BondNum);
                    
                    Target1 = [5 epitope];
                    Target2 = [epitope 5];
                 
                    for i=1:length(Bondsidx)
                        bonds = data.Bonds{Bondsidx(i)};
                        idx5 = find( (bonds(:,4)==5) | (bonds(:,5)==5) ); %% All attachments of Ab to epitope
                        if(length(idx5) & isempty(Kon))
                            MFET = dt*data.timestep(Bondsidx(i));
                            Kon =  MFET;
                            Kon2ndarm = NaN;
                        end
                        
                        if(length(idx5) & (LatestMFPT==0)) % The lastest arrival time of one arm to the epitope
                            LatestMFPT = dt*data.timestep(Bondsidx(i));
                        end
                        
                        Arm1idx = 0;
                        Arm2idx = 0;
                        SpikeArm1 = 0;
                        SpikeArm2 = 0;

                        OneArmBindTime = [OneArmBindTime ; data.timestep(Bondsidx(i)) ];
                        for b=1:length(idx5)
                            bond = bonds(idx5(b),:);
                            
                            if( isequal(bond(4:5),Target1) | isequal(bond(4:5),Target2) );
                                BondEnergy = [BondEnergy , bond(7)];
                                BondDist = [BondDist , bond(6)];

                                Position = find(bond(4:5)==5);
                                EpitopePosition = find(bond(4:5)==epitope);
                                idxArm = bond(Position+1);
                                idxSpike = bond(EpitopePosition+1);
                                if(Arm1idx==0)
                                    Arm1idx = idxArm;
                                    OneArmcount = OneArmcount+1;
                                    SpikeArm1 = idxSpike;
                                elseif(idxArm~=Arm1idx)
                                    
                                    if(isnan(Kon2ndarm))
                                        
                                        % find the latest time 1 Ab was
                                        % bound
                                        
                                        d1ArmBound = diff(OneArmBindTime);
                                        DTSim = min(d1ArmBound);
                                        dtTime = 1.05*min(d1ArmBound);
                                        d1ArmArrival = find(d1ArmBound>dtTime);
                                        if(length(d1ArmArrival))
                                            LatestArrvialtime = OneArmBindTime(d1ArmArrival(end)+1);
                                        else
                                            LatestArrvialtime = MFET;
                                        end

                                        MFET2arm = dt*(data.timestep(Bondsidx(i))-LatestArrvialtime);
                                        if(MFET2arm==0)
                                            MFET2arm = dt*DTSim;
                                            Kon2ndarm = MFET2arm;
                                        else
                                            if(Prev1ArmMFPT~=LatestArrvialtime)
                                                if(isnan(Kon2ndarm))
                                                    Kon2ndarm = MFET2arm;
                                                else
                                                    Kon2ndarm = [Kon2ndarm ; MFET2arm];
                                                end
                                            end

                                        end
                                        Prev1ArmMFPT = LatestArrvialtime;
                                    end

                                    TwoArmcount = TwoArmcount+1;
                                    
                                    
                                    
                                end
                                
                                Boundscount = Boundscount+1;

                            end
                        end
                    end
                end
                KonN_epitope = [KonN_epitope ; Kon];
                KonN_2ndArm_epitope = [KonN_2ndArm_epitope ; Kon2ndarm];
                
            end
            end

            KonArr(n,ep) = mean(KonN_epitope);
            idxnotnan = find(~isnan(KonN_2ndArm_epitope));
            Kon2ndArmArr(n,ep) = mean(KonN_2ndArm_epitope(idxnotnan));
            
            KonArrCell{n,ep} = KonN_epitope;
            Kon2ndArmArrCell{n,ep} = KonN_2ndArm_epitope;
            
            BondEnergytArr(n,ep) = mean(BondEnergy);
            BondDistArr(n,ep) = mean(BondDist);
            TwoArmDistBoundArr(n,ep) = mean(TwoArmDistBound);
            TwoArmDistBoundSameSpikeArr(n,ep) = TwoArmDistBoundSameSpike;

            BoundscountArr(n,ep) = Boundscount/BoundsTotcountArr(n,ep);
            OneArmcountArr(n,ep) = OneArmcount/BoundsTotcountArr(n,ep);
            TwoArmcountArr(n,ep) = TwoArmcount/BoundsTotcountArr(n,ep);
            
            
            
        
        
    end
end

P1Arm = zeros(length(NArr),length(epitopeArr));
P2Arm = zeros(length(NArr),length(epitopeArr));

for i=1:length(NArr)
    for j=1:length(epitopeArr)
        OneArmMFPT = KonArrCell{i,j};
        idxnotnan = find(~isnan(OneArmMFPT));
        P1Arm(i,j) = length(idxnotnan)/BoundscountArrNorm(i,j);
        
        TwoArmMFPT = Kon2ndArmArrCell{i,j};
        idxnotnan = find(~isnan(TwoArmMFPT));
        P2Arm(i,j) = length(idxnotnan)/BoundscountArrNorm(i,j);
        
    end
end

kOn1Arm = zeros(length(NArr),length(epitopeArr));
kOn2Arm = zeros(length(NArr),length(epitopeArr));


for i=1:length(NArr)
    for j=1:length(epitopeArr)
        OneArmMFPT = KonArrCell{i,j};
        idxnotnan = find(~isnan(OneArmMFPT));
        if(length(idxnotnan))
            kOn1Arm(i,j) = 1/mean(OneArmMFPT(idxnotnan));
        end
        
        TwoArmMFPT = Kon2ndArmArrCell{i,j};
        idxnotnan = find(~isnan(TwoArmMFPT));
        if(length(idxnotnan))
            kOn2Arm(i,j) = 1/mean(TwoArmMFPT(idxnotnan));
        end
    end
end

flag_str = ['virusha40Spikes'];

kOn1Arm_HA = kOn1Arm;
kOn2Arm_HA = kOn2Arm;
OneArmcountArr_HA = OneArmcountArr;
TwoArmcountArr_HA = TwoArmcountArr;
P1Arm_HA = P1Arm;
P2Arm_HA = P2Arm;

save([flag_str,'.mat'],'NArr','BoundscountArr','OneArmcountArr','TwoArmcountArr','KonArr','Kon2ndArmArr','KonArrCell',...
    'Kon2ndArmArrCell','BoundscountArrNorm','BoundsTotcountArr','BondEnergytArr','BondDistArr','TwoArmDistBoundArr','TwoArmDistBoundSameSpikeArr',...
    'P1Arm','P2Arm','kOn1Arm','kOn2Arm','P1Arm_HA','P2Arm_HA','kOn1Arm_HA','kOn2Arm_HA','OneArmcountArr_HA','TwoArmcountArr_HA','-v7.3');


