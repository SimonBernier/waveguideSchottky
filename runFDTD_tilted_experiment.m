function runFDTD_tilted_experiment(Array_id)
%This accuracy requires about 2300MB of RAM and about 12hours to run for
%the longest runTimes (v=0.73c, a=0.5mm, b=6mm). It should prevent
%numerical instabilities at plasma frequencies of 3THz. The runtime is
%adaptive and is calculated as T=t0+(b-a)/v+(7e-3-b)/cSi. This experiment
%uses a Gaussian profile in z, increases zskip to 15 such that we can
%resolve the detector position within dz=6.8 microns. This corresponds to a
%FFT calculation with a max frequency of 6.40THz. The case v=c is made to
%have tskip=3 instead of 6 in an effort to resolve the pulse while still
%minimizing RAM and physical storage.

maxf = 6e12;
v_list = [0.73 0.86 1.0 1.07 1.49];
fp_list = [0.5 1 3 7]*1e12;
Profile_list = 1:13;
[v,fp,Profile]=meshgrid(v_list,fp_list,Profile_list);
v=reshape(v,1,[]); fp=reshape(fp,1,[]); Profile=reshape(Profile,1,[]);
Es=1e4; Lz = 8e-3;

% length(v) = 260 in this case
if v(Array_id)==v_list(3)
    %run special code for v=c case
    tic
    FDTD_tilted(maxf, fp(Array_id), v(Array_id), Es, Lz, Profile(Array_id), 3, 15)
    fprintf('run %d : %0.1f s\n', Array_id, toc)
else %run normal code
    tic
    FDTD_tilted(maxf, fp(Array_id), v(Array_id), Es, Lz, Profile(Array_id), 6, 15)
    fprintf('run %d : %0.1f s\n', Array_id, toc)
end