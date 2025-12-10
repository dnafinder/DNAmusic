function DNAmusic3(varargin)
%DNAMUSIC3 - The life music
%Audification (or sonification) is the technique of using the sense of hearing
%to analyse data. The advantage of audification over visualisation in data
%analysis is that sound has the property that when different notes are played
%togheter they can still be individually heard: in vision colours blend to form
%new colours. DNA and proteins maps naturally onto musical sequences. Several
%algorithms were proposed to translate DNA and proteins into music. This
%function uses the algorithm proposed by Rie Takahashi and Jeffrey Miller
%(Genome Biology 2007; 8(5):405).
%
% Syntax: 	DNAmusic3(mRNAid,wn)
%
% Inputs:
%   mRNAid - this is the id of the messanger that you want to translate
%            deposited on the NCBI database.
%            This is a string.
%            (default = 'NM_005218').
%   wn     - this is the duration (in second) of a whole note (default = 2).
%
% Example:
%   DNAmusic3
%
% Of course longer is the mRNA and longer (very much longer) will be the time to
% convert it into music and much higher memory needed...
%
% Created by Giuseppe Cardillo
% giuseppe.cardillo.75@gmail.com
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). DNAmusic3 - The life music by Takahashi and Miller
% https://github.com/dnafinder/DNAmusic

% -------------------------
% Input error handling
% -------------------------
args = cell(varargin);
nu = numel(args);
if nu > 2
    error('DNAmusic3 accepts max 2 input arguments')
end

default.values = {'NM_005218'; 2};
default.values(1:nu) = args;
[mRNAid, semibreve] = deal(default.values{:});

if nu >= 1 && ~(ischar(mRNAid) || (isstring(mRNAid) && isscalar(mRNAid)))
    error('DNAmusic3 requires a string as mRNAid')
end
mRNAid = char(mRNAid);

if nu == 2 && ( ~isscalar(semibreve) || ~isfinite(semibreve) || ~isnumeric(semibreve) )
    error('DNAmusic3 requires a scalar, numeric and finite WN value.')
end

% -------------------------
% Retrieve informations
% -------------------------
try
    S = genbankread([mRNAid '.GBK']);
catch
    disp(['Query the NCBI database to retrieve the ' mRNAid ' mRNA'])
    try
        S = getgenbank(mRNAid);
        disp('Data retrieved')
    catch ME
        error('DNAmusic3:NCBIRetrievalFailed', ...
            'Unable to retrieve GenBank record for %s. %s', mRNAid, ME.message);
    end
end

% Defensive checks
if ~isfield(S,'Sequence') || ~isfield(S,'CDS') || ...
        ~isfield(S.CDS,'indices') || ~isfield(S.CDS,'translation')
    error('DNAmusic3:MalformedGenBank', ...
        'GenBank record for %s does not contain required CDS information.', mRNAid);
end

% -------------------------
% Keep only the useful informations
% -------------------------
mRNAsequence = upper(S.Sequence);   % mRNA sequence
mRNAcoding  = S.CDS.indices;        % the coding portion
protein     = S.CDS.translation;    % the protein sequence

if isfield(S,'Definition')
    disp(['Sonification process for ' S.Definition])
else
    disp('Sonification process for retrieved mRNA.')
end

if numel(mRNAcoding) < 2 || any(~isfinite(mRNAcoding))
    error('DNAmusic3:InvalidCDS', 'Invalid CDS indices in GenBank record for %s.', mRNAid);
end

startCDS = mRNAcoding(1);
endCDS   = mRNAcoding(2);

if startCDS < 1 || endCDS > numel(mRNAsequence) || startCDS >= endCDS
    error('DNAmusic3:InvalidCDSRange', 'CDS indices out of range for %s.', mRNAid);
end

% -------------------------
% Audio parameters
% -------------------------
fs = 8192;          % historical sample rate used by the original implementation
fadeTime = 0.01;    % sec (anti-click fade)

% -------------------------
% Codon frequencies are used to set duration of the notes
% -------------------------
CDS = double(mRNAsequence(startCDS:endCDS)); % ASCII codes for CDS bases

if mod(numel(CDS),3) ~= 0
    CDS = CDS(1:floor(numel(CDS)/3)*3); % defensive truncation
end

% Create an Nx3 matrix: each row is a codon (ASCII)
cd = reshape(CDS, 3, numel(CDS)/3)';

% Embedded Standard DNA code matrix
code = getDNAcode();

% Find codons in the matrix
[tf, codonIdx] = ismember(cd, code(:,1:3), 'rows');
if any(~tf)
    error('DNAmusic3:UnknownCodon', ...
        'One or more CDS codons were not found in the embedded DNA code table.');
end

% Pick the codon frequency (column 5)
cf = code(codonIdx, 5)';

% The CDS includes a stop codon; the protein does not.
lp = length(protein);

% Defensive alignment: if cf has one extra element, remove the last one.
if numel(cf) == lp + 1
    cf(end) = [];
elseif numel(cf) ~= lp
    % Keep behavior safe: align to the shorter length
    n = min(numel(cf), lp);
    cf = cf(1:n);
    protein = protein(1:n);
    lp = n;
end

% Assign note durations by codon frequency classes (as in original implementation)
d = zeros(1, lp);
L = zeros(1, lp);

idx = (cf > 0) & (cf < 11);
d(idx) = semibreve/8;  L(idx) = length(0:1/fs:semibreve/8);

idx = (cf >= 11) & (cf < 21);
d(idx) = semibreve/4;  L(idx) = length(0:1/fs:semibreve/4);

idx = (cf >= 21) & (cf < 30);
d(idx) = semibreve/2;  L(idx) = length(0:1/fs:semibreve/2);

idx = (cf >= 30);
d(idx) = semibreve;    L(idx) = length(0:1/fs:semibreve);

% Any zero-length entries fallback to shortest class
idx = (L == 0);
if any(idx)
    d(idx) = semibreve/8;
    L(idx) = length(0:1/fs:semibreve/8);
end

% -------------------------
% Convert the protein into triads using logical indexing
% -------------------------
disp('Protein sonification')
fprintf('\tConvert %i amino acids into a notes array\n', lp)

Proteinmusic = cell(lp, 3); % each AA mapped to a triad

Proteinmusic(protein=='W',1) = {'C3'};  Proteinmusic(protein=='W',2) = {'E3'};  Proteinmusic(protein=='W',3) = {'G3'};
Proteinmusic(protein=='M',1) = {'D3'};  Proteinmusic(protein=='M',2) = {'F#3'}; Proteinmusic(protein=='M',3) = {'A3'};
Proteinmusic(protein=='P',1) = {'E3'};  Proteinmusic(protein=='P',2) = {'G#3'}; Proteinmusic(protein=='P',3) = {'B3'};
Proteinmusic(protein=='H',1) = {'F3'};  Proteinmusic(protein=='H',2) = {'A3'};  Proteinmusic(protein=='H',3) = {'C4'};
Proteinmusic(protein=='Y',1) = {'G3'};  Proteinmusic(protein=='Y',2) = {'B3'};  Proteinmusic(protein=='Y',3) = {'D4'};
Proteinmusic(protein=='F',1) = {'B3'};  Proteinmusic(protein=='F',2) = {'D4'};  Proteinmusic(protein=='F',3) = {'G4'};
Proteinmusic(protein=='L',1) = {'A3'};  Proteinmusic(protein=='L',2) = {'C#4'}; Proteinmusic(protein=='L',3) = {'E4'};
Proteinmusic(protein=='I',1) = {'C#4'}; Proteinmusic(protein=='I',2) = {'E4'};  Proteinmusic(protein=='I',3) = {'A4'};
Proteinmusic(protein=='V',1) = {'B3'};  Proteinmusic(protein=='V',2) = {'D#4'}; Proteinmusic(protein=='V',3) = {'F#4'};
Proteinmusic(protein=='A',1) = {'D#4'}; Proteinmusic(protein=='A',2) = {'F#4'}; Proteinmusic(protein=='A',3) = {'B4'};
Proteinmusic(protein=='C',1) = {'C4'};  Proteinmusic(protein=='C',2) = {'E4'};  Proteinmusic(protein=='C',3) = {'G4'};
Proteinmusic(protein=='G',1) = {'D4'};  Proteinmusic(protein=='G',2) = {'F#4'}; Proteinmusic(protein=='G',3) = {'A4'};
Proteinmusic(protein=='T',1) = {'E4'};  Proteinmusic(protein=='T',2) = {'G#4'}; Proteinmusic(protein=='T',3) = {'B4'};
Proteinmusic(protein=='S',1) = {'G#4'}; Proteinmusic(protein=='S',2) = {'B4'};  Proteinmusic(protein=='S',3) = {'E5'};
Proteinmusic(protein=='Q',1) = {'F4'};  Proteinmusic(protein=='Q',2) = {'A4'};  Proteinmusic(protein=='Q',3) = {'C5'};
Proteinmusic(protein=='N',1) = {'A4'};  Proteinmusic(protein=='N',2) = {'C5'};  Proteinmusic(protein=='N',3) = {'F5'};
Proteinmusic(protein=='E',1) = {'G4'};  Proteinmusic(protein=='E',2) = {'B4'};  Proteinmusic(protein=='E',3) = {'D5'};
Proteinmusic(protein=='D',1) = {'B4'};  Proteinmusic(protein=='D',2) = {'D5'};  Proteinmusic(protein=='D',3) = {'G5'};
Proteinmusic(protein=='R',1) = {'A4'};  Proteinmusic(protein=='R',2) = {'C#5'}; Proteinmusic(protein=='R',3) = {'E5'};
Proteinmusic(protein=='K',1) = {'C#5'}; Proteinmusic(protein=='K',2) = {'E5'};  Proteinmusic(protein=='K',3) = {'A5'};

% Defensive fill for any unassigned cells
emptyCell = cellfun(@isempty, Proteinmusic);
if any(emptyCell(:))
    Proteinmusic(emptyCell) = {'r'};
end

fprintf('\tDone\n')
fprintf('\tConvert notes array into frequencies array\n')

% Build sample index ranges per amino acid
idxstop  = cumsum(L);
idxstart = [1, idxstop(1:end-1)+1];

% Convert into tunes
fprintf('\t\tConvert 1st root\n')
first = zeros(idxstop(end), 1);
for k = 1:lp
    first(idxstart(k):idxstop(k)) = fnote(Proteinmusic{k,1}, d(k), fs, fadeTime);
end
fprintf('\t\tDone\n')

fprintf('\t\tConvert 3rd major\n')
third = zeros(idxstop(end), 1);
for k = 1:lp
    third(idxstart(k):idxstop(k)) = fnote(Proteinmusic{k,2}, d(k), fs, fadeTime);
end
fprintf('\t\tDone\n')

fprintf('\t\tConvert 5th perfect\n')
fifth = zeros(idxstop(end), 1);
for k = 1:lp
    % Correctly use the 3rd column for the fifth note of the triad
    fifth(idxstart(k):idxstop(k)) = fnote(Proteinmusic{k,3}, d(k), fs, fadeTime);
end
fprintf('\t\tDone\n')

disp('Mix and play')
disp('to stop press ctrl-c')
disp(' ')

mix = first + third + fifth;
playWithProgress(mix, fs, "Protein sonification (triads)")

end

% =========================================================================
% Local functions
% =========================================================================
function y = fnote(str, dur, fs, fadeTime)
%FNOTE Convert note string into waveform for the given duration.
% str: e.g., 'A3', 'F#4', or 'r' (rest)

t = 0:1/fs:dur;

if isempty(str) || strcmp(str,'r')
    y = zeros(size(t));
    return
end

% Piano key frequencies
a = ['A';' ';'B';'C';' ';'D';' ';'E';'F';' ';'G']; % keyboard
iv = find(a == str(1), 1, 'first');
if isempty(iv)
    y = zeros(size(t));
    return
end

oct = str2double(str(end)) * 12;
alt = (length(str) == 3);
key = iv + alt + oct;
f = 440 * (2^(1/12))^(key - 49);

y = sin((2*pi*f) .* t);

% Anti-click fade
nFade = floor(fs * fadeTime) + 1;
nFade = min(nFade, length(y));

if nFade > 1
    fadein  = linspace(0, 1, nFade);
    fadeout = linspace(1, 0, nFade);
    y(1:nFade) = y(1:nFade) .* fadein;
    y(end-nFade+1:end) = y(end-nFade+1:end) .* fadeout;
end
end

function playWithProgress(x, fs, sectionTitle)
%PLAYWITHPROGRESS Play audio with a simple progress bar using audioplayer.

if isempty(x)
    return
end

x = x(:);
p = audioplayer(x, fs);
totalSamples = numel(x);

wb = waitbar(0, sectionTitle, 'Name', 'DNAmusic progress');

p.TimerPeriod = 0.1;
p.TimerFcn = @(~,~) updateWB();
p.StopFcn  = @(~,~) closeWB();

    function updateWB()
        if isvalid(p) && ishandle(wb)
            frac = min(1, max(0, p.CurrentSample / totalSamples));
            waitbar(frac, wb, sprintf('%s (%.0f%%)', sectionTitle, frac*100));
        end
    end

    function closeWB()
        if ishandle(wb)
            close(wb);
        end
    end

try
    playblocking(p);
catch ME
    if ishandle(wb)
        close(wb);
    end
    rethrow(ME);
end
end

function code = getDNAcode()
%GETDNACODE Embedded codon table used by DNAmusic3.
%
% Columns:
% 1-3: DNA codon in ASCII (A=65, C=67, G=71, T=84)
% 4  : Amino acid in ASCII (includes stop)
% 5  : Codon frequency (used to set note duration)
% 6  : Additional normalized/relative value (kept for fidelity)

code = [ ...
65 65 65 75 24.4 0.519563; ...
65 65 67 78 19.1 1; ...
65 65 71 75 31.9 0.685819; ...
65 65 84 78 17 0.607579; ...
65 67 65 84 15.1 0.183405; ...
65 67 67 84 18.9 0.220049; ...
65 67 71 84 6.1 0.242054; ...
65 67 84 84 13.1 0.305623; ...
65 71 65 82 12.2 0.183374; ...
65 71 67 83 19.5 0.244499; ...
65 71 71 82 12 0.211491; ...
65 71 84 83 12.1 0.144254; ...
65 84 65 73 7.5 0.152855; ...
65 84 67 73 20.8 0.552567; ...
65 84 71 77 22 0.611247; ...
65 84 84 73 16 0.572127; ...
67 65 65 81 12.3 0.336186; ...
67 65 67 72 15.1 0.336186; ...
67 65 71 81 34.2 0.749389; ...
67 65 84 72 10.9 0.19835; ...
67 67 65 80 16.9 0.213967; ...
67 67 67 80 19.8 0.220049; ...
67 67 71 80 6.9 0.190709; ...
67 67 84 80 17.5 0.305623; ...
67 71 65 82 6.2 0.183395; ...
67 71 67 82 10.4 0.154034; ...
67 71 71 82 11.4 0.211491; ...
67 71 84 82 4.5 0.213936; ...
67 84 65 76 7.2 0.091724; ...
67 84 67 76 19.6 0.264059; ...
67 84 71 76 39.6 0.334963; ...
67 84 84 76 13.2 0.366748; ...
71 65 65 69 29 0.397311; ...
71 65 67 68 25.1 0.580685; ...
71 65 71 69 39.6 0.52445; ...
71 65 84 68 21.8 0.342604; ...
71 67 65 65 15.8 0.27515; ...
71 67 67 65 27.7 0.638142; ...
71 67 71 65 7.4 0.240831; ...
71 67 84 65 18.4 0.886308; ...
71 71 65 71 16.5 0.275061; ...
71 71 67 71 22.2 0.458435; ...
71 71 71 71 16.5 0.301956; ...
71 71 84 71 10.8 0.270477; ...
71 84 65 86 7.1 0.152845; ...
71 84 67 86 14.5 0.242054; ...
71 84 71 86 28.1 0.537897; ...
71 84 84 86 11 0.336186; ...
84 65 65 88 1 0; ...
84 65 67 89 15.3 0.449878; ...
84 65 71 88 0.8 0; ...
84 65 84 89 12.2 0.283007; ...
84 67 65 83 12.2 0.152845; ...
84 67 67 83 17.7 0.242054; ...
84 67 71 83 4.4 0.171149; ...
84 67 84 83 15.2 0.336186; ...
84 71 65 88 1.6 0; ...
84 71 67 67 12.6 0.91687; ...
84 71 71 87 13.2 0.275061; ...
84 71 84 67 10.6 0.540954; ...
84 84 65 76 7.7 0.213936; ...
84 84 67 70 20.3 0.366748; ...
84 84 71 76 12.9 0.282396; ...
84 84 84 70 17.6 0.216381 ...
];
end
