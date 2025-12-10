function DNAmusic2(varargin)
%DNAMUSIC2 - The life music
%Audification (or sonification) is the technique of using the sense of hearing
%to analyse data. The advantage of audification over visualisation in data
%analysis is that sound has the property that when different notes are played
%togheter they can still be individually heard: in vision colours blend to form
%new colours. DNA and proteins maps naturally onto musical sequences. Several
%algorithms were proposed to translate DNA and proteins into music. This
%function use the algorithm proposed by Ross D. King and Colin G. Angus
%(PM - Protein Music - Cabios applications notes 1996; 12(3):251-252)
%
% Syntax: 	DNAmusic2(mRNAid,wn)
%
% Inputs:
%   mRNAid - this is the id of the messanger that you want to translate
%            deposited on the NCBI database.
%            This is a string.
%            (default = 'NM_005218').
%   wn     - this is the duration (in second) of a whole note (default = 1).
%
% Example:
%   DNAmusic2
%
% Of course longer is the mRNA and longer (very much longer) will be the time to
% convert it into music and much higher memory needed...
%
% Created by Giuseppe Cardillo
% giuseppe.cardillo.75@gmail.com
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). DNAmusic2 - The life music by King and Angus
% https://github.com/dnafinder/DNAmusic

% -------------------------
% Input error handling
% -------------------------
args = cell(varargin);
nu = numel(args);
if nu > 2
    error('DNAmusic2 accepts max 2 input arguments')
end

default.values = {'NM_005218'; 1};
default.values(1:nu) = args;
[mRNAid, semibreve] = deal(default.values{:});

if nu >= 1 && ~(ischar(mRNAid) || (isstring(mRNAid) && isscalar(mRNAid)))
    error('DNAmusic2 requires a string as mRNAid')
end
mRNAid = char(mRNAid);

if nu == 2 && ( ~isscalar(semibreve) || ~isfinite(semibreve) || ~isnumeric(semibreve) )
    error('DNAmusic2 requires a scalar, numeric and finite WN value.')
end

% -------------------------
% Retrieve informations
% -------------------------
try
    % Try to see if a genbank file was saved on disk
    S = genbankread([mRNAid '.GBK']);
catch
    % Else retrieve informations from the NCBI database
    disp(['Query the NCBI database to retrieve the ' mRNAid ' mRNA'])
    try
        S = getgenbank(mRNAid);
        disp('Data retrieved')
    catch ME
        error('DNAmusic2:NCBIRetrievalFailed', ...
            'Unable to retrieve GenBank record for %s. %s', mRNAid, ME.message);
    end
end

% Defensive checks
if ~isfield(S,'Sequence') || ~isfield(S,'CDS') || ...
        ~isfield(S.CDS,'indices') || ~isfield(S.CDS,'translation')
    error('DNAmusic2:MalformedGenBank', ...
        'GenBank record for %s does not contain required CDS information.', mRNAid);
end

% -------------------------
% Keep only the useful informations
% -------------------------
mRNAsequence = upper(S.Sequence);     % mRNA sequence
mRNAcoding  = S.CDS.indices;          % the coding portion
protein     = S.CDS.translation;      % the protein sequence

if isfield(S,'Definition')
    disp(['Sonification process for ' S.Definition])
else
    disp('Sonification process for retrieved mRNA.')
end

% -------------------------
% Audio parameters
% -------------------------
fs = 8192;          % historical sample rate used by the original implementation
fadeTime = 0.01;    % sec (anti-click fade)

% -------------------------
% Start the mRNA sonification process
% 1) Convert all nucleotides into notes according to:
% PM - Protein Music - Cabios applications notes 1996; 12(3):251-252
% -------------------------
dur = semibreve/4;
t = 0:1/fs:dur;
L = length(t);

disp('mRNA sonification')
fprintf('....Convert %i nucleotides into a notes array\n', length(mRNAsequence))

mRNAmusic = cell(length(mRNAsequence),1); % array preallocation
mRNAmusic(mRNAsequence=='A') = {'A4'};
mRNAmusic(mRNAsequence=='C') = {'C3'};
mRNAmusic(mRNAsequence=='G') = {'G3'};
mRNAmusic(mRNAsequence=='T') = {'E3'};

% Split the music in three parts:
if numel(mRNAcoding) < 2 || any(~isfinite(mRNAcoding))
    error('DNAmusic2:InvalidCDS', 'Invalid CDS indices in GenBank record for %s.', mRNAid);
end

startCDS = mRNAcoding(1);
endCDS   = mRNAcoding(2);

if startCDS < 1 || endCDS > numel(mRNAmusic) || startCDS >= endCDS
    error('DNAmusic2:InvalidCDSRange', 'CDS indices out of range for %s.', mRNAid);
end

% Prelude (5' UTR)
mRNAmusicprelude = mRNAmusic(1:startCDS-1);
% Main theme (CDS)
mRNAmusicmaintheme = mRNAmusic(startCDS:endCDS);
% Finale (3' UTR)
mRNAmusicfinale = mRNAmusic(endCDS+1:end);

disp('....Done')
disp('....Convert notes array into frequencies array')

% Convert into tunes (mRNA)
% Prelude
l = length(mRNAmusicprelude);
lhprelude = zeros(L*l,1);
for k = 1:l
    lhprelude((k-1)*L+1:k*L) = fnote(mRNAmusicprelude{k}, t, fs, fadeTime);
end

% Main theme
l = length(mRNAmusicmaintheme);
lhmaintheme = zeros(L*l,1);
for k = 1:l
    lhmaintheme((k-1)*L+1:k*L) = fnote(mRNAmusicmaintheme{k}, t, fs, fadeTime);
end

% Finale
l = length(mRNAmusicfinale);
lhfinale = zeros(L*l,1);
for k = 1:l
    lhfinale((k-1)*L+1:k*L) = fnote(mRNAmusicfinale{k}, t, fs, fadeTime);
end

disp('....Done')

% -------------------------
% Protein sonification
% -------------------------
dur = semibreve/8;
t = 0:1/fs:dur;
L = length(t);

disp('Protein sonification')
fprintf('....Convert %i amino acids into a notes array\n', length(protein))

% The array must be 1 cell longer because the last codon is a stop codon:
% so there is not amino acid (insert a rest).
Proteinmusic = cell(length(protein)+1,6);

% Default: leave as empty then fill; initialize proline as rests
idxP = (protein == 'P');
Proteinmusic(idxP,:) = {'r'};

% I, V, L
idx = find(protein=='I' | protein=='V' | protein=='L');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'C1' 'C1' 'r' 'r' 'G1' 'G1'};
end

% W, Y
idx = find(protein=='W' | protein=='Y');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'D1' 'D1' 'C1' 'C1' 'A2' 'A2'};
end

% M, C, A, G
idx = find(protein=='M' | protein=='C' | protein=='A' | protein=='G');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'C1' 'C1' 'r' 'r' 'r' 'r'};
end

% S, Q, N
idx = find(protein=='S' | protein=='Q' | protein=='N');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'A2' 'A2' 'r' 'r' 'r' 'r'};
end

% T, E
idx = find(protein=='T' | protein=='E');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'A2' 'A2' 'r' 'r' 'C1' 'C1'};
end

% D
idx = find(protein=='D');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'A2' 'A2' 'r' 'r' 'F1' 'F1'};
end

% H
idx = find(protein=='H');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'E1' 'E1' 'A2' 'C1' 'D1' 'F1'};
end

% K
idx = find(protein=='K');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'A2' 'A2' 'C1' 'C1' 'F1' 'E1'};
end

% R
idx = find(protein=='R');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'A2' 'A2' 'E1' 'E1' 'F1' 'F1'};
end

% F
idx = find(protein=='F');
for k = 1:numel(idx)
    Proteinmusic(idx(k),:) = {'D1' 'D1' 'r' 'r' 'C1' 'C1'};
end

% Fill any still-empty rows with rests (defensive)
emptyRow = cellfun(@isempty, Proteinmusic(:,1));
Proteinmusic(emptyRow,:) = {'r'};

% Append stop codon rests
Proteinmusic(end,:) = {'r'};

Proteinmusic = Proteinmusic';
Proteinmusic = Proteinmusic(:);

disp('....Done')
disp('....Convert notes array into frequencies array')

% Convert into tunes (protein main theme)
l = length(Proteinmusic);
rhmaintheme = zeros(L*l,1);
for k = 1:l
    rhmaintheme((k-1)*L+1:k*L) = fnote(Proteinmusic{k}, t, fs, fadeTime);
end

disp('Done')

% -------------------------
% Mix and play with progress
% -------------------------
disp('Mix and play')
disp('to stop press ctrl-c')
disp(' ')

disp('Prelude (mRNA 5'' UTR)')
playWithProgress(lhprelude, fs, "Prelude (mRNA 5' UTR)")

disp('Main theme (CDS and protein)')
u = min(length(rhmaintheme), length(lhmaintheme));
mixMain = lhmaintheme(1:u) + rhmaintheme(1:u);
playWithProgress(mixMain, fs, "Main theme (CDS and protein)")

disp('Finale (mRNA 3'' UTR)')
playWithProgress(lhfinale, fs, "Finale (mRNA 3' UTR)")

end

% =========================================================================
% Local functions
% =========================================================================
function y = fnote(str, t, fs, fadeTime)
%FNOTE Convert note string into waveform for the given time vector.
% str: e.g., 'A3', 'F#4', or 'r' (rest)

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
