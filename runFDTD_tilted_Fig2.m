function runFDTD_tilted_Fig2(ArrayId)

maxf = 6e12;
v = 0.86;
Profile = 1:7;
Es=1e4; Lz = 8e-3;

tic
FDTD_tilted_Fig2(maxf, v, Es, Lz, Profile(ArrayId), 6, 12)
fprintf('run %d : %0.1f s\n', Array_id, toc)
