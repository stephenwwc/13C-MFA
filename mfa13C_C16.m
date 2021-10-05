function obj = mfa13C_C16(var)   
%% precusor information
% this version delete v18 reverse AcCoA -> PYR
 r1 = var(1);
 r2 = var(2);
 r3 = var(3);
 r4 = var(4);
 r5 = var(5);
 r6 = var(6);
 r7 = var(7);
 r8 = var(8);
 r9 = var(9);
 r10 = var(10);
 r11 = var(11);
 r12 = var(12);
 r13 = var(13);
 r14 = var(14);
 r15 = var(15);
 r16 = var(16);
 r17 = var(17);
 r18 = var(18);
 br1 = var(19);

% precusor demand ratio 
% here we can modify the community
GramPos = [0.47 5.01 5.53 6.28 3.27 4.28 2.54]';
GramNeg = [0.32 9.29 11.63 13.2 5.01 7.44 4.89]';
Fungi = [0.18 1.30 1.45 1.72 1.04 1.08 0.6]';
% if use gramneg
% Biomassf = br1*GramNeg; % if choose GramNeg
% Biomassf = br1*GramPos; % if choose GramPos
% Biomassf = br1*Fungi; % if choose Fungi

% % if use community, 
Bplus = 4.5;
Bminus = 4.5;
Fungus = 1;
Community = [Bplus Bminus Fungus]/(Bplus+Bminus+Fungus).*[GramPos, GramNeg, Fungi];
Community = sum(Community, 2);
Biomassf = br1*Community;  % attention, the commuinity chosen

% obtaine other bimass fluxes br2 to br8
br2 = Biomassf(1);
br3 = Biomassf(2);
br4 = Biomassf(3);
br5 = Biomassf(4);
br6 = Biomassf(5);
br7 = Biomassf(6);
br8 = Biomassf(7);



%% define the atom mapping matrix 
% R1, 'glucose -> glucose-P'
MatrixGLU_G6P = diag([1 1 1 1 1 1]);
% R2 AND R13, 'G6P -> F6P
MatrixG6P_F6P = diag([1 1 1 1 1 1]);
MatrixF6P_G6P = diag([1 1 1 1 1 1]);

% R3 AND R14, F6P -> GAP(abc) + GAP(def)
MatrixF6P_GAPinr3 = 1/2*[flip(diag([1 1 1])), diag([1 1 1])];
% R14, reverse reaction GAP(abc) + GAP(ABC) -> F6P(cbaABC)
MatrixGAP_F6Pinr14a = [flip(diag([1 1 1])); zeros(3)];
MatrixGAP_F6Pinr14b = [zeros(3,3); diag([1 1 1])];

% R4 AND R15,'GAP -> PYR + 22NADH + 2 ATP
MatrixGAP_PYR = diag([1 1 1]);
MatrixPYR_GAP = diag([1 1 1]);

% R5 PYR -> AcCoA + co2 + NADPH
MatrixPYR_AcCoA = [0 1 0; 0 0 1];
MatrixPYR_CO2 = [1 0 0];

% R6 AcCoA + OAA -> ICIT
MatrixAcCoA_ICIT = [1 0;0 1; zeros(4,2)];
MatrixOAA_ICIT = [zeros(2,4);0 1 0 0;0 0 1 0;0 0 0 1; 1 0 0 0];

% R7 ICIT > KG + CO2 + NADH
MatrixICIT_AKG = [diag([1 1 1 1 1]), zeros(5,1)];
MatrixICIT_CO2 = [0 0 0 0 0 1];

% R8 aKG > 1/2 OAA  + 1/2 OAA + CO2 + NADH + ATP + FADH2  !!! pay attetion
MatrixAKG_OAA1 = [diag([1 1 1 1]), zeros(4,1)];
MatrixAKG_OAA2 = [flip(diag([1 1 1 1])), zeros(4,1)];
MatrixAKG_OAA = [MatrixAKG_OAA1 + MatrixAKG_OAA2]/2;
MatrixAKG_CO2 = [0 0 0 0 1];

% R9 G6P  > RU5P + CO2 + 2 NADPH
 MatrixG6P_RU5P = [[0 0 0 0 0]', diag([1 1 1 1 1])];
 MatrixG6P_CO2 = [1 0 0 0 0 0];

% R10 RU5P + RU5P -> S7P + GAP
MatrixRU5P_S7P1 = [zeros(2,5); diag([1 1 1 1 1])];
MatrixRU5P_S7P2 = [diag([1 1 0 0 0]); zeros(2,5)];
MatrixRU5P_GAPinr10 = [zeros(3,2), diag([1 1 1])];

% R11 S7P + GAP -> F6P + E4P;
MatrixS7P_E4P = [zeros(4,3), diag([1 1 1 1])];
MatrixS7P_F6Pinr11 = [diag([1 1 1 0 0 0]), zeros(6,1)];
MatrixGAP_F6Pinr11 = [zeros(3,3); diag([1 1 1])];

% R12 E4P + RU5P -> F6P + GAP
MatrixE4P_F6Pinr12 = [zeros(2,4); diag([1 1 1 1])];
MatrixRU5P_F6Pinr12 = [diag([1 1 0 0 0]); zeros(1, 5)];
MatrixRU5P_GAPinr12 = [zeros(3,2), diag([1 1 1])];

% R16 ABD r17 PYR + CO2 + ATP -> OAA
MatrixPYR_OAA = [diag([1 1 1]); zeros(1,3)];
MatrixCO2_OAA = [0 0 0 1]';
% r17
MatrixOAA_PYR = MatrixPYR_OAA';
MatrixOAA_CO2 = MatrixCO2_OAA';

% R18 G6P(abcdef) -> PYR(abc) + GAP(def)
MatrixG6P_PYRinr18 = [diag([1 1 1]),zeros(3)];
MatrixG6P_GAPinr18= [zeros(3), diag([1 1 1])];

%% 1-13C-GLU %%
GLU1 = [1 0 0 0 0 0]';
% set the intial value
GLU1 = [1 0 0 0 0 0]';
G6P1 = [1 0 0 0 0 0]';
F6P1 = [0.01 0 0 0 0 0]';
GAP1  = [0.01 0 0]';
PYR1 = [0.01 0 0]';
AcCoA1 = [0.01 0]';
ICIT1 = [0.01 0 0 0 0 0]';
AKG1 = [0.01 0 0 0 0]';
OAA1 = [0.01 0 0 0]';
RU5P1 = [0.01 0 0 0 0]';
S7P1 = [0.01 0 0 0 0 0 0]';
E4P1 = [0.01 0 0 0]';

% calculate the each fragment MID
iter = 50;
for k = 2:iter
G6P1(:,k) = (MatrixGLU_G6P*GLU1 * r1 + MatrixF6P_G6P * F6P1(:,k-1) * r13 ) / (r2 + r9 + r18 + br1); % add new r18

F6P1(:,k) = 1/(r3 + r13 + br2)*(MatrixG6P_F6P * G6P1(:,k) * r2 + (MatrixGAP_F6Pinr14a * GAP1(:,k-1)...
    + MatrixGAP_F6Pinr14b * GAP1(:,k-1))*r14...
    + (MatrixS7P_F6Pinr11*S7P1(:,k-1) + MatrixGAP_F6Pinr11 * GAP1(:,k-1)) * r11...
    + (MatrixRU5P_F6Pinr12*RU5P1(:,k-1) + MatrixE4P_F6Pinr12 * E4P1(:,k-1))* r12);

GAP1(:,k) = 1/(r11 + r4 + r14 + br3)* (MatrixF6P_GAPinr3 * F6P1(:,k) * 2*r3...
    + MatrixPYR_GAP * PYR1(:,k-1) * r15...
    + MatrixRU5P_GAPinr12 * RU5P1(:,k-1) * r12...
    + MatrixRU5P_GAPinr10 * RU5P1(:,k-1) * r10...
    + MatrixG6P_GAPinr18* G6P1(:,k-1)* r18); % add new r18

PYR1(:,k) = 1/(r5 + r15 + r16 + br4)...
     *(MatrixGAP_PYR * GAP1(:,k) * r4...
    + MatrixOAA_PYR * OAA1(:,k-1) * r17...
    + MatrixG6P_PYRinr18 * G6P1(:,k-1) *r18); % add new r18

AcCoA1(:,k) = 1/(r6 + br5)*(MatrixPYR_AcCoA * PYR1(:,k) * r5);

ICIT1(:,k) = 1/r7 * (MatrixAcCoA_ICIT * AcCoA1(:,k) + MatrixOAA_ICIT * OAA1(:,k-1)) * r6;

AKG1(:,k) =1/(r8 + br6)* MatrixICIT_AKG * ICIT1(:,k) * r7;

OAA1(:,k) = 1/(r6 + r17 + br7)*(MatrixAKG_OAA * AKG1(:,k) * r8 + MatrixPYR_OAA * PYR1(:,k) * r16);

RU5P1(:,k) = 1/(2*r10 + r12 + br8)* MatrixG6P_RU5P * G6P1(:,k) * r9;

S7P1(:,k) = 1/r11* (MatrixRU5P_S7P1 * RU5P1(:,k) + MatrixRU5P_S7P2 * RU5P1(:,k))* r10;

E4P1(:,k) = 1/r12* (MatrixS7P_E4P * S7P1(:,k))* r11;

G6P1error = abs(G6P1(:,k)-G6P1(:,k-1));
F6P1error = abs(F6P1(:,k)-F6P1(:,k-1));
GAP1error = abs(GAP1(:,k)-GAP1(:,k-1));
PYR1error = abs(PYR1(:,k)-PYR1(:,k-1));
AcCoA1error = abs(AcCoA1(:,k)-AcCoA1(:,k-1));
ICIT1error = abs(ICIT1(:,k)-ICIT1(:,k-1));
AKG1error = abs(AKG1(:,k)-AKG1(:,k-1));
OAA1error = abs(OAA1(:,k)-OAA1(:,k-1));
RU5P1error = abs(RU5P1(:,k)-RU5P1(:,k-1));
S7P1error = abs(S7P1(:,k)-S7P1(:,k-1));
E4P1error = abs(E4P1(:,k)-E4P1(:,k-1));
% CO21error = abs(CO21(:,k)-CO21(:,k-1));


end
totalerror1 = sum(G6P1error) + sum(F6P1error) + sum(F6P1error) + sum(GAP1error)...
    +sum(PYR1error) + sum(AcCoA1error) + sum(ICIT1error) + sum(AKG1error)...
    +sum(OAA1error) + sum(RU5P1error) + sum(S7P1error) + sum(E4P1error);

%% 2-13C-GLU
GLU2 = [0 1 0 0 0 0]';

% set the intial value
G6P2 = [1 0 0 0 0 0]';
F6P2 = [0.01 0 0 0 0 0]';
GAP2  = [0.01 0 0]';
PYR2 = [0.01 0 0]';
AcCoA2 = [0.01 0]';
ICIT2 = [0.01 0 0 0 0 0]';
AKG2 = [0.01 0 0 0 0]';
OAA2 = [0.01 0 0 0]';
RU5P2 = [0.01 0 0 0 0]';
S7P2 = [0.01 0 0 0 0 0 0]';
E4P2 = [0.01 0 0 0]';




%% calculate the each fragment MID
iter = 50;
for k = 2:iter
G6P2(:,k) = (MatrixGLU_G6P*GLU2 * r1 + MatrixF6P_G6P * F6P2(:,k-1) * r13 ) / (r2 + r9 + r18 + br1); % add new r18

F6P2(:,k) = 1/(r3 + r13 + br2)*(MatrixG6P_F6P * G6P2(:,k) * r2 + (MatrixGAP_F6Pinr14a * GAP2(:,k-1)...
    + MatrixGAP_F6Pinr14b * GAP2(:,k-1))*r14...
    + (MatrixS7P_F6Pinr11*S7P2(:,k-1) + MatrixGAP_F6Pinr11 * GAP2(:,k-1)) * r11...
    + (MatrixRU5P_F6Pinr12*RU5P2(:,k-1) + MatrixE4P_F6Pinr12 * E4P2(:,k-1))* r12);

GAP2(:,k) = 1/(r11 + r4 + r14 + br3)* (MatrixF6P_GAPinr3 * F6P2(:,k) * 2*r3...
    + MatrixPYR_GAP * PYR2(:,k-1) * r15...
    + MatrixRU5P_GAPinr12 * RU5P2(:,k-1) * r12...
    + MatrixRU5P_GAPinr10 * RU5P2(:,k-1) * r10...
    + MatrixG6P_GAPinr18* G6P2(:,k-1)* r18); % add new r18

PYR2(:,k) = 1/(r5 + r15 + r16 + br4)...
     *(MatrixGAP_PYR * GAP2(:,k) * r4...
    + MatrixOAA_PYR * OAA2(:,k-1) * r17...
    + MatrixG6P_PYRinr18 * G6P2(:,k-1) *r18); % add new r18

AcCoA2(:,k) = 1/(r6 + br5)*(MatrixPYR_AcCoA * PYR2(:,k) * r5);

ICIT2(:,k) = 1/r7 * (MatrixAcCoA_ICIT * AcCoA2(:,k) + MatrixOAA_ICIT * OAA2(:,k-1)) * r6;

AKG2(:,k) =1/(r8 + br6)* MatrixICIT_AKG * ICIT2(:,k) * r7;

OAA2(:,k) = 1/(r6 + r17 + br7)*(MatrixAKG_OAA * AKG2(:,k) * r8 + MatrixPYR_OAA * PYR2(:,k) * r16);

RU5P2(:,k) = 1/(2*r10 + r12 + br8)* MatrixG6P_RU5P * G6P2(:,k) * r9;

S7P2(:,k) = 1/r11* (MatrixRU5P_S7P1 * RU5P2(:,k) + MatrixRU5P_S7P2 * RU5P2(:,k))* r10;

E4P2(:,k) = 1/r12* (MatrixS7P_E4P * S7P2(:,k))* r11;


G6P2error = abs(G6P2(:,k)-G6P2(:,k-1));
F6P2error = abs(F6P2(:,k)-F6P2(:,k-1));
GAP2error = abs(GAP2(:,k)-GAP2(:,k-1));
PYR2error = abs(PYR2(:,k)-PYR2(:,k-1));
AcCoA2error = abs(AcCoA2(:,k)-AcCoA2(:,k-1));
ICIT2error = abs(ICIT2(:,k)-ICIT2(:,k-1));
AKG2error = abs(AKG2(:,k)-AKG2(:,k-1));
OAA2error = abs(OAA2(:,k)-OAA2(:,k-1));
RU5P2error = abs(RU5P2(:,k)-RU5P2(:,k-1));
S7P2error = abs(S7P2(:,k)-S7P2(:,k-1));
E4P2error = abs(E4P2(:,k)-E4P2(:,k-1));



end
totalerror2 = sum(G6P2error) + sum(F6P2error) + sum(F6P2error) + sum(GAP2error)...
    +sum(PYR2error) + sum(AcCoA2error) + sum(ICIT2error) + sum(AKG2error)...
    +sum(OAA2error) + sum(RU5P2error) + sum(S7P2error) + sum(E4P2error);

%% 3-13C-GLU
GLU3 = [0 0 1 0 0 0]';

% set the intial value
G6P3 = [1 0 0 0 0 0]';
F6P3 = [0.01 0 0 0 0 0]';
GAP3  = [0.01 0 0]';
PYR3 = [0.01 0 0]';
AcCoA3 = [0.01 0]';
ICIT3 = [0.01 0 0 0 0 0]';
AKG3 = [0.01 0 0 0 0]';
OAA3 = [0.01 0 0 0]';
RU5P3 = [0.01 0 0 0 0]';
S7P3 = [0.01 0 0 0 0 0 0]';
E4P3 = [0.01 0 0 0]';
% CO23 = [0.01]';

% calculate the each fragment MID
iter = 50;
for k = 2:iter
G6P3(:,k) = (MatrixGLU_G6P*GLU3 * r1 + MatrixF6P_G6P * F6P3(:,k-1) * r13 ) / (r2 + r9 + r18 + br1); % add new r18

F6P3(:,k) = 1/(r3 + r13 + br2)*(MatrixG6P_F6P * G6P3(:,k) * r2 + (MatrixGAP_F6Pinr14a * GAP3(:,k-1)...
    + MatrixGAP_F6Pinr14b * GAP3(:,k-1))*r14...
    + (MatrixS7P_F6Pinr11*S7P3(:,k-1) + MatrixGAP_F6Pinr11 * GAP3(:,k-1)) * r11...
    + (MatrixRU5P_F6Pinr12*RU5P3(:,k-1) + MatrixE4P_F6Pinr12 * E4P3(:,k-1))* r12);

GAP3(:,k) = 1/(r11 + r4 + r14 + br3)* (MatrixF6P_GAPinr3 * F6P3(:,k) * 2*r3...
    + MatrixPYR_GAP * PYR3(:,k-1) * r15...
    + MatrixRU5P_GAPinr12 * RU5P3(:,k-1) * r12...
    + MatrixRU5P_GAPinr10 * RU5P3(:,k-1) * r10...
    + MatrixG6P_GAPinr18* G6P3(:,k-1)* r18); % add new r18

PYR3(:,k) = 1/(r5 + r15 + r16 + br4)...
     *(MatrixGAP_PYR * GAP3(:,k) * r4...
    + MatrixOAA_PYR * OAA3(:,k-1) * r17...
    + MatrixG6P_PYRinr18 * G6P3(:,k-1) *r18); % add new r18

AcCoA3(:,k) = 1/(r6 + br5)*(MatrixPYR_AcCoA * PYR3(:,k) * r5);

ICIT3(:,k) = 1/r7 * (MatrixAcCoA_ICIT * AcCoA3(:,k) + MatrixOAA_ICIT * OAA3(:,k-1)) * r6;

AKG3(:,k) =1/(r8 + br6)* MatrixICIT_AKG * ICIT3(:,k) * r7;

OAA3(:,k) = 1/(r6 + r17 + br7)*(MatrixAKG_OAA * AKG3(:,k) * r8 + MatrixPYR_OAA * PYR3(:,k) * r16);

RU5P3(:,k) = 1/(2*r10 + r12 + br8)* MatrixG6P_RU5P * G6P3(:,k) * r9;

S7P3(:,k) = 1/r11* (MatrixRU5P_S7P1 * RU5P3(:,k) + MatrixRU5P_S7P2 * RU5P3(:,k))* r10;

E4P3(:,k) = 1/r12* (MatrixS7P_E4P * S7P3(:,k))* r11;
% 
% CO23(:,k) = 1/(1+r16)* (MatrixG6P_CO2 * G6P3(:,k) * r9...
%     + MatrixPYR_CO2 * PYR3(:,k) * r5...
%     + MatrixICIT_CO2 * ICIT3(:,k) * r7...
%     + MatrixAKG_CO2 * AKG3(:,k) * r8...
%     + MatrixOAA_CO2 * OAA3(:,k) * r17);

G6P3error = abs(G6P3(:,k)-G6P3(:,k-1));
F6P3error = abs(F6P3(:,k)-F6P3(:,k-1));
GAP3error = abs(GAP3(:,k)-GAP3(:,k-1));
PYR3error = abs(PYR3(:,k)-PYR3(:,k-1));
AcCoA3error = abs(AcCoA3(:,k)-AcCoA3(:,k-1));
ICIT3error = abs(ICIT3(:,k)-ICIT3(:,k-1));
AKG3error = abs(AKG3(:,k)-AKG3(:,k-1));
OAA3error = abs(OAA3(:,k)-OAA3(:,k-1));
RU5P3error = abs(RU5P3(:,k)-RU5P3(:,k-1));
S7P3error = abs(S7P3(:,k)-S7P3(:,k-1));
E4P3error = abs(E4P3(:,k)-E4P3(:,k-1));
% CO23error = abs(CO23(:,k)-CO23(:,k-1));


end
totalerror3 = sum(G6P3error) + sum(F6P3error) + sum(F6P3error) + sum(GAP3error)...
    +sum(PYR3error) + sum(AcCoA3error) + sum(ICIT3error) + sum(AKG3error)...
    +sum(OAA3error) + sum(RU5P3error) + sum(S7P3error) + sum(E4P3error);


%%  4-13C-glu stimulation
GLU4 = [0 0 0 1 0 0]'; %%

% set the intial value
G6P4 = [0 0 0 1 0 0]';
F6P4 = [0.01 0 0 0 0 0]';
GAP4  = [0.01 0 0]';
PYR4 = [0.01 0 0]';
AcCoA4 = [0.01 0]';
ICIT4 = [0.01 0 0 0 0 0]';
AKG4 = [0.01 0 0 0 0]';
OAA4 = [0.01 0 0 0]';
RU5P4 = [0.01 0 0 0 0]';
S7P4 = [0.01 0 0 0 0 0 0]';
E4P4 = [0.01 0 0 0]';
% CO24 = [0.01]';



% calculate the each fragment MID
iter = 50;
for k = 2:iter
G6P4(:,k) = (MatrixGLU_G6P*GLU4 * r1 + MatrixF6P_G6P * F6P4(:,k-1) * r13 ) / (r2 + r9 + r18 + br1); % add new r18

F6P4(:,k) = 1/(r3 + r13 + br2)*(MatrixG6P_F6P * G6P4(:,k) * r2 + (MatrixGAP_F6Pinr14a * GAP4(:,k-1)...
    + MatrixGAP_F6Pinr14b * GAP4(:,k-1))*r14...
    + (MatrixS7P_F6Pinr11*S7P4(:,k-1) + MatrixGAP_F6Pinr11 * GAP4(:,k-1)) * r11...
    + (MatrixRU5P_F6Pinr12*RU5P4(:,k-1) + MatrixE4P_F6Pinr12 * E4P4(:,k-1))* r12);

GAP4(:,k) = 1/(r11 + r4 + r14 + br3)* (MatrixF6P_GAPinr3 * F6P4(:,k) * 2*r3...
    + MatrixPYR_GAP * PYR4(:,k-1) * r15...
    + MatrixRU5P_GAPinr12 * RU5P4(:,k-1) * r12...
    + MatrixRU5P_GAPinr10 * RU5P4(:,k-1) * r10...
    + MatrixG6P_GAPinr18* G6P4(:,k-1)* r18); % add new r18

PYR4(:,k) = 1/(r5 + r15 + r16 + br4)...
     *(MatrixGAP_PYR * GAP4(:,k) * r4...
    + MatrixOAA_PYR * OAA4(:,k-1) * r17...
    + MatrixG6P_PYRinr18 * G6P4(:,k-1) *r18); % add new r18

AcCoA4(:,k) = 1/(r6 + br5)*(MatrixPYR_AcCoA * PYR4(:,k) * r5);

ICIT4(:,k) = 1/r7 * (MatrixAcCoA_ICIT * AcCoA4(:,k) + MatrixOAA_ICIT * OAA4(:,k-1)) * r6;

AKG4(:,k) =1/(r8 + br6)* MatrixICIT_AKG * ICIT4(:,k) * r7;

OAA4(:,k) = 1/(r6 + r17 + br7)*(MatrixAKG_OAA * AKG4(:,k) * r8 + MatrixPYR_OAA * PYR4(:,k) * r16);

RU5P4(:,k) = 1/(2*r10 + r12 + br8)* MatrixG6P_RU5P * G6P4(:,k) * r9;

S7P4(:,k) = 1/r11* (MatrixRU5P_S7P1 * RU5P4(:,k) + MatrixRU5P_S7P2 * RU5P4(:,k))* r10;

E4P4(:,k) = 1/r12* (MatrixS7P_E4P * S7P4(:,k))* r11;
% CO24(:,k) = 1/(1+r16)* (MatrixG6P_CO2 * G6P4(:,k) * r9...
%     + MatrixPYR_CO2 * PYR4(:,k) * r5...
%     + MatrixICIT_CO2 * ICIT4(:,k) * r7...
%     + MatrixAKG_CO2 * AKG4(:,k) * r8...
%     + MatrixOAA_CO2 * OAA4(:,k) * r17);

G6P4error = abs(G6P4(:,k)-G6P4(:,k-1));
F6P4error = abs(F6P4(:,k)-F6P4(:,k-1));
GAP4error = abs(GAP4(:,k)-GAP4(:,k-1));
PYR4error = abs(PYR4(:,k)-PYR4(:,k-1));
AcCoA4error = abs(AcCoA4(:,k)-AcCoA4(:,k-1));
ICIT4error = abs(ICIT4(:,k)-ICIT4(:,k-1));
AKG4error = abs(AKG4(:,k)-AKG4(:,k-1));
OAA4error = abs(OAA4(:,k)-OAA4(:,k-1));
RU5P4error = abs(RU5P4(:,k)-RU5P4(:,k-1));
S7P4error = abs(S7P4(:,k)-S7P4(:,k-1));
E4P4error = abs(E4P4(:,k)-E4P4(:,k-1));



end
totalerror4 = sum(G6P4error) + sum(F6P4error) + sum(F6P4error) + sum(GAP4error)...
    +sum(PYR4error) + sum(AcCoA4error) + sum(ICIT4error) + sum(AKG4error)...
    +sum(OAA4error) + sum(RU5P4error) + sum(S7P4error) + sum(E4P4error);

%%  5-13C-glucose stimulation
GLU5 = [0 0 0 0 1 0]'; %%

% set the intial value
G6P5 = [1 0 0 0 0 0]';
F6P5 = [0.01 0 0 0 0 0]';
GAP5  = [0.01 0 0]';
PYR5 = [0.01 0 0]';
AcCoA5 = [0.01 0]';
ICIT5 = [0.01 0 0 0 0 0]';
AKG5 = [0.01 0 0 0 0]';
OAA5 = [0.01 0 0 0]';
RU5P5 = [0.01 0 0 0 0]';
S7P5 = [0.01 0 0 0 0 0 0]';
E4P5 = [0.01 0 0 0]';

% calculate the each fragment MID
iter = 50;
for k = 2:iter
G6P5(:,k) = (MatrixGLU_G6P*GLU5 * r1 + MatrixF6P_G6P * F6P5(:,k-1) * r13 ) / (r2 + r9 + r18 + br1); % add new r18

F6P5(:,k) = 1/(r3 + r13 + br2)*(MatrixG6P_F6P * G6P5(:,k) * r2 + (MatrixGAP_F6Pinr14a * GAP5(:,k-1)...
    + MatrixGAP_F6Pinr14b * GAP5(:,k-1))*r14...
    + (MatrixS7P_F6Pinr11*S7P5(:,k-1) + MatrixGAP_F6Pinr11 * GAP5(:,k-1)) * r11...
    + (MatrixRU5P_F6Pinr12*RU5P5(:,k-1) + MatrixE4P_F6Pinr12 * E4P5(:,k-1))* r12);

GAP5(:,k) = 1/(r11 + r4 + r14 + br3)* (MatrixF6P_GAPinr3 * F6P5(:,k) * 2*r3...
    + MatrixPYR_GAP * PYR5(:,k-1) * r15...
    + MatrixRU5P_GAPinr12 * RU5P5(:,k-1) * r12...
    + MatrixRU5P_GAPinr10 * RU5P5(:,k-1) * r10...
    + MatrixG6P_GAPinr18* G6P5(:,k-1)* r18); % add new r18

PYR5(:,k) = 1/(r5 + r15 + r16 + br4)...
     *(MatrixGAP_PYR * GAP5(:,k) * r4...
    + MatrixOAA_PYR * OAA5(:,k-1) * r17...
    + MatrixG6P_PYRinr18 * G6P5(:,k-1) *r18); % add new r18

AcCoA5(:,k) = 1/(r6 + br5)*(MatrixPYR_AcCoA * PYR5(:,k) * r5);

ICIT5(:,k) = 1/r7 * (MatrixAcCoA_ICIT * AcCoA5(:,k) + MatrixOAA_ICIT * OAA5(:,k-1)) * r6;

AKG5(:,k) =1/(r8 + br6)* MatrixICIT_AKG * ICIT5(:,k) * r7;

OAA5(:,k) = 1/(r6 + r17 + br7)*(MatrixAKG_OAA * AKG5(:,k) * r8 + MatrixPYR_OAA * PYR5(:,k) * r16);

RU5P5(:,k) = 1/(2*r10 + r12 + br8)* MatrixG6P_RU5P * G6P5(:,k) * r9;

S7P5(:,k) = 1/r11* (MatrixRU5P_S7P1 * RU5P5(:,k) + MatrixRU5P_S7P2 * RU5P5(:,k))* r10;

E4P5(:,k) = 1/r12* (MatrixS7P_E4P * S7P5(:,k))* r11;

G6P5error = abs(G6P5(:,k)-G6P5(:,k-1));
F6P5error = abs(F6P5(:,k)-F6P5(:,k-1));
GAP5error = abs(GAP5(:,k)-GAP5(:,k-1));
PYR5error = abs(PYR5(:,k)-PYR5(:,k-1));
AcCoA5error = abs(AcCoA5(:,k)-AcCoA5(:,k-1));
ICIT5error = abs(ICIT5(:,k)-ICIT5(:,k-1));
AKG5error = abs(AKG5(:,k)-AKG5(:,k-1));
OAA5error = abs(OAA5(:,k)-OAA5(:,k-1));
RU5P5error = abs(RU5P5(:,k)-RU5P5(:,k-1));
S7P5error = abs(S7P5(:,k)-S7P5(:,k-1));
E4P5error = abs(E4P5(:,k)-E4P5(:,k-1));

end
totalerror5 = sum(G6P5error) + sum(F6P5error) + sum(F6P5error) + sum(GAP5error)...
    +sum(PYR5error) + sum(AcCoA5error) + sum(ICIT5error) + sum(AKG5error)...
    +sum(OAA5error) + sum(RU5P5error) + sum(S7P5error) + sum(E4P5error);

%%  6-13C-glucose stimulation

GLU6 = [0 0 0 0 0 1]'; %%

% set the intial value
G6P6 = [1 0 0 0 0 0]';
F6P6 = [0.01 0 0 0 0 0]';
GAP6  = [0.01 0 0]';
PYR6 = [0.01 0 0]';
AcCoA6 = [0.01 0]';
ICIT6 = [0.01 0 0 0 0 0]';
AKG6 = [0.01 0 0 0 0]';
OAA6 = [0.01 0 0 0]';
RU5P6 = [0.01 0 0 0 0]';
S7P6 = [0.01 0 0 0 0 0 0]';
E4P6 = [0.01 0 0 0]';

% calculate the each fragment MID
iter = 50;
for k = 2:iter
G6P6(:,k) = (MatrixGLU_G6P*GLU6 * r1 + MatrixF6P_G6P * F6P6(:,k-1) * r13 ) / (r2 + r9 + r18 + br1); % add new r18

F6P6(:,k) = 1/(r3 + r13 + br2)*(MatrixG6P_F6P * G6P6(:,k) * r2 + (MatrixGAP_F6Pinr14a * GAP6(:,k-1)...
    + MatrixGAP_F6Pinr14b * GAP6(:,k-1))*r14...
    + (MatrixS7P_F6Pinr11*S7P6(:,k-1) + MatrixGAP_F6Pinr11 * GAP6(:,k-1)) * r11...
    + (MatrixRU5P_F6Pinr12*RU5P6(:,k-1) + MatrixE4P_F6Pinr12 * E4P6(:,k-1))* r12);

GAP6(:,k) = 1/(r11 + r4 + r14 + br3)* (MatrixF6P_GAPinr3 * F6P6(:,k) * 2*r3...
    + MatrixPYR_GAP * PYR6(:,k-1) * r15...
    + MatrixRU5P_GAPinr12 * RU5P6(:,k-1) * r12...
    + MatrixRU5P_GAPinr10 * RU5P6(:,k-1) * r10...
    + MatrixG6P_GAPinr18* G6P6(:,k-1)* r18); % add new r18

PYR6(:,k) = 1/(r5 + r15 + r16 + br4)...
     *(MatrixGAP_PYR * GAP6(:,k) * r4...
    + MatrixOAA_PYR * OAA6(:,k-1) * r17...
    + MatrixG6P_PYRinr18 * G6P6(:,k-1) *r18); % add new r18

AcCoA6(:,k) = 1/(r6 + br5)*(MatrixPYR_AcCoA * PYR6(:,k) * r5);

ICIT6(:,k) = 1/r7 * (MatrixAcCoA_ICIT * AcCoA6(:,k) + MatrixOAA_ICIT * OAA6(:,k-1)) * r6;

AKG6(:,k) =1/(r8 + br6)* MatrixICIT_AKG * ICIT6(:,k) * r7;

OAA6(:,k) = 1/(r6 + r17 + br7)*(MatrixAKG_OAA * AKG6(:,k) * r8 + MatrixPYR_OAA * PYR6(:,k) * r16);

RU5P6(:,k) = 1/(2*r10 + r12 + br8)* MatrixG6P_RU5P * G6P6(:,k) * r9;

S7P6(:,k) = 1/r11* (MatrixRU5P_S7P1 * RU5P6(:,k) + MatrixRU5P_S7P2 * RU5P6(:,k))* r10;

E4P6(:,k) = 1/r12* (MatrixS7P_E4P * S7P6(:,k))* r11;

G6P6error = abs(G6P6(:,k)-G6P6(:,k-1));
F6P6error = abs(F6P6(:,k)-F6P6(:,k-1));
GAP6error = abs(GAP6(:,k)-GAP6(:,k-1));
PYR6error = abs(PYR6(:,k)-PYR6(:,k-1));
AcCoA6error = abs(AcCoA6(:,k)-AcCoA6(:,k-1));
ICIT6error = abs(ICIT6(:,k)-ICIT6(:,k-1));
AKG6error = abs(AKG6(:,k)-AKG6(:,k-1));
OAA6error = abs(OAA6(:,k)-OAA6(:,k-1));
RU5P6error = abs(RU5P6(:,k)-RU5P6(:,k-1));
S7P6error = abs(S7P6(:,k)-S7P6(:,k-1));
E4P6error = abs(E4P6(:,k)-E4P6(:,k-1));



end
totalerror6 = sum(G6P6error) + sum(F6P6error) + sum(F6P6error) + sum(GAP6error)...
    +sum(PYR6error) + sum(AcCoA6error) + sum(ICIT6error) + sum(AKG6error)...
    +sum(OAA6error) + sum(RU5P6error) + sum(S7P6error) + sum(E4P6error);


%%  U-13C6-GLUCOSE
GLUU = [1 1 1 1 1 1]'; %%

% set the intial value
G6PU = [1 1 1 1 1 1]';
F6PU = [0.01 0 0 0 0 0]';
GAPU  = [0.01 0 0]';
PYRU = [0.01 0 0]';
AcCoAU = [0.01 0]';
ICITU = [0.01 0 0 0 0 0]';
AKGU = [0.01 0 0 0 0]';
OAAU = [0.01 0 0 0]';
RU5PU = [0.01 0 0 0 0]';
S7PU = [0.01 0 0 0 0 0 0]';
E4PU = [0.01 0 0 0]';

% calculate the each fragment MID
iter = 50;
for k = 2:iter
G6PU(:,k) = (MatrixGLU_G6P*GLUU * r1 + MatrixF6P_G6P * F6PU(:,k-1) * r13 ) / (r2 + r9 + r18 + br1); % add new r18

F6PU(:,k) = 1/(r3 + r13 + br2)*(MatrixG6P_F6P * G6PU(:,k) * r2 + (MatrixGAP_F6Pinr14a * GAPU(:,k-1)...
    + MatrixGAP_F6Pinr14b * GAPU(:,k-1))*r14...
    + (MatrixS7P_F6Pinr11*S7PU(:,k-1) + MatrixGAP_F6Pinr11 * GAPU(:,k-1)) * r11...
    + (MatrixRU5P_F6Pinr12*RU5PU(:,k-1) + MatrixE4P_F6Pinr12 * E4PU(:,k-1))* r12);

GAPU(:,k) = 1/(r11 + r4 + r14 + br3)* (MatrixF6P_GAPinr3 * F6PU(:,k) * 2*r3...
    + MatrixPYR_GAP * PYRU(:,k-1) * r15...
    + MatrixRU5P_GAPinr12 * RU5PU(:,k-1) * r12...
    + MatrixRU5P_GAPinr10 * RU5PU(:,k-1) * r10...
    + MatrixG6P_GAPinr18* G6PU(:,k-1)* r18); % add r18

PYRU(:,k) = 1/(r5 + r15 + r16 + br4)...
     *(MatrixGAP_PYR * GAPU(:,k) * r4...
    + MatrixOAA_PYR * OAAU(:,k-1) * r17...
    + MatrixG6P_PYRinr18 * G6PU(:,k-1) *r18); % add new r18

AcCoAU(:,k) = 1/(r6 + br5)*(MatrixPYR_AcCoA * PYRU(:,k) * r5);

ICITU(:,k) = 1/r7 * (MatrixAcCoA_ICIT * AcCoAU(:,k) + MatrixOAA_ICIT * OAAU(:,k-1)) * r6;

AKGU(:,k) =1/(r8 + br6)* MatrixICIT_AKG * ICITU(:,k) * r7;

OAAU(:,k) = 1/(r6 + r17 + br7)*(MatrixAKG_OAA * AKGU(:,k) * r8 + MatrixPYR_OAA * PYRU(:,k) * r16);

RU5PU(:,k) = 1/(2*r10 + r12 + br8)* MatrixG6P_RU5P * G6PU(:,k) * r9;

S7PU(:,k) = 1/r11* (MatrixRU5P_S7P1 * RU5PU(:,k) + MatrixRU5P_S7P2 * RU5PU(:,k))* r10;

E4PU(:,k) = 1/r12* (MatrixS7P_E4P * S7PU(:,k))* r11;

G6PUerror = abs(G6PU(:,k)-G6PU(:,k-1));
F6PUerror = abs(F6PU(:,k)-F6PU(:,k-1));
GAPUerror = abs(GAPU(:,k)-GAPU(:,k-1));
PYRUerror = abs(PYRU(:,k)-PYRU(:,k-1));
AcCoAUerror = abs(AcCoAU(:,k)-AcCoAU(:,k-1));
ICITUerror = abs(ICITU(:,k)-ICITU(:,k-1));
AKGUerror = abs(AKGU(:,k)-AKGU(:,k-1));
OAAUerror = abs(OAAU(:,k)-OAAU(:,k-1));
RU5PUerror = abs(RU5PU(:,k)-RU5PU(:,k-1));
S7PUerror = abs(S7PU(:,k)-S7PU(:,k-1));
E4PUerror = abs(E4PU(:,k)-E4PU(:,k-1));



end
totalerrorU = sum(G6PUerror) + sum(F6PUerror) + sum(F6PUerror) + sum(GAPUerror)...
    +sum(PYRUerror) + sum(AcCoAUerror) + sum(ICITUerror) + sum(AKGUerror)...
    +sum(OAAUerror) + sum(RU5PUerror) + sum(S7PUerror) + sum(E4PUerror);



INPUT = [0.069	0.185	0.147	0.038	0.228	0.236	0.052	0.219	0.129	0.033	0.283	0.176]; % corrected with Uniformaly 13C-alebled


% for PF
% propyl 
Propyl1 = [AcCoA1(:,50); AcCoA1(1,50)];
Propyl2 = [AcCoA2(:,50); AcCoA2(1,50)];
Propyl3 = [AcCoA3(:,50); AcCoA3(1,50)];
Propyl4 = [AcCoA4(:,50); AcCoA4(1,50)];
Propyl5 = [AcCoA5(:,50); AcCoA5(1,50)];
Propyl6 = [AcCoA6(:,50); AcCoA6(1,50)];
PropylU = [AcCoAU(:,50); AcCoAU(1,50)];


% %!!! we should input the stimulated FL_AcCoA
MZ74 = sum([AcCoA1(:,50) AcCoA2(:,50) AcCoA3(:,50) AcCoA4(:,50) AcCoA5(:,50) AcCoA6(:,50)])/sum(AcCoAU(:,50));
MZ87 = sum([Propyl1 Propyl2 Propyl3 Propyl4 Propyl5 Propyl6])/sum(PropylU);
SIMU = [MZ74 MZ87];
% % necessary to output
obj = sum((SIMU - INPUT).^2);

