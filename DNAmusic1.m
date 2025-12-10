function DNAmusic1(varargin)
%DNAMUSIC1 - The life music
%Audification (or sonification) is the technique of using the sense of hearing
%to analyse data. The advantage of audification over visualisation in data
%analysis is that sound has the property that when different notes are played
%togheter they can still be individually heard: in vision colours blend to form
%new colours. DNA and proteins maps naturally onto musical sequences. Several
%algorithms were proposed to translate DNA and proteins into music. This
%function use the algorithm proposed by Nobuo Munakata
%(http://www.toshima.ne.jp/~edogiku/)
%
% Syntax: 	DNAmusic1(mRNAid,wn)
%
% Inputs:
%   mRNAid - this is the id of the messanger that you want to translate
%            deposited on the NCBI database
%            (http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene).
%            This is a string.
%            (default = 'NM_005218').
%   wn     - this is the duration (in second) of a whole note (default = 2).
%
% Example:
%   DNAmusic1
%
% Of course longer is the mRNA and longer (very much longer) will be the time to
% convert it into music and much higher memory needed...
%
% Created by Giuseppe Cardillo
% giuseppe.cardillo.75@gmail.com
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). DNAmusic1 - The life music by Nobuo Munakata
% https://github.com/dnafinder/DNAmusic

% -------------------------
% Input error handling
% -------------------------
args = cell(varargin);
nu = numel(args);
if nu > 2
    error('DNAmusic1 accepts max 2 input arguments')
end

default.values = {'NM_005218'; 2};
default.values(1:nu) = args;
[mRNAid, semibreve] = deal(default.values{:});

if nu >= 1 && ~ischar(mRNAid)
    error('DNAmusic1 requires a string as mRNAid')
end

if nu == 2 && ( ~isscalar(semibreve) || ~isfinite(semibreve) || ~isnumeric(semibreve) )
    error('DNAmusic1 requires a scalar, numeric and finite WN value.')
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
        error('DNAmusic1:NCBIRetrievalFailed', ...
            'Unable to retrieve GenBank record for %s. %s', mRNAid, ME.message);
    end
end

% Basic field checks (defensive)
if ~isfield(S,'Sequence') || ~isfield(S,'CDS') || ~isfield(S.CDS,'indices') || ~isfield(S.CDS,'translation')
    error('DNAmusic1:MalformedGenBank', ...
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
fs = 8192;                 % historical sample rate used by the original implementation
fadeTime = 0.01;           % sec (anti-click fade)
dur = semibreve/8;         % duration of each note (as in original DNAmusic1)

t = 0:1/fs:dur;
L = length(t);             % length of each note vector

% -------------------------
% Start the mRNA sonification process
% 1) Convert all nucleotides into notes according to Munakata mapping
% -------------------------
disp('mRNA sonification')
fprintf('....Convert %i nucleotides into a notes array\n', length(mRNAsequence))

mRNAmusic = cell(length(mRNAsequence),1); % array preallocation
mRNAmusic(mRNAsequence=='A') = {'A3'};
mRNAmusic(mRNAsequence=='C') = {'E3'};
mRNAmusic(mRNAsequence=='G') = {'D3'};
mRNAmusic(mRNAsequence=='T') = {'G3'};

% -------------------------
% 2) Split the music in three parts:
% -------------------------
if numel(mRNAcoding) < 2 || any(~isfinite(mRNAcoding))
    error('DNAmusic1:InvalidCDS', 'Invalid CDS indices in GenBank record for %s.', mRNAid);
end

% Prelude (5' UTR)
startCDS = mRNAcoding(1);
endCDS   = mRNAcoding(2);

if startCDS < 1 || endCDS > numel(mRNAmusic) || startCDS >= endCDS
    error('DNAmusic1:InvalidCDSRange', 'CDS indices out of range for %s.', mRNAid);
end

mRNAmusicprelude = mRNAmusic(1:startCDS-1);

% Main theme (CDS)
mRNAmusicmaintheme = mRNAmusic(startCDS:endCDS);

% Ensure divisibility by 3 before reshape (defensive)
if mod(length(mRNAmusicmaintheme),3) ~= 0
    % Keep behavior safe: truncate trailing bases if malformed record
    mRNAmusicmaintheme = mRNAmusicmaintheme(1:floor(length(mRNAmusicmaintheme)/3)*3);
end

mRNAmusicmaintheme = reshape(mRNAmusicmaintheme, 3, length(mRNAmusicmaintheme)/3);
mRNAmusicmaintheme(4,:) = {'r'}; % insert rest to articulate codon group
mRNAmusicmaintheme = mRNAmusicmaintheme(:);

% Finale (3' UTR)
mRNAmusicfinale = mRNAmusic(endCDS+1:end);

disp('....Done')
disp('....Convert notes array into frequencies array')

% -------------------------
% Convert into tunes
% -------------------------
% Prelude
l = length(mRNAmusicprelude);
lhprelude = zeros(L*l,1); % array preallocation
for k = 1:l
    lhprelude((k-1)*L+1:k*L) = fnote(mRNAmusicprelude{k}, t, fs, fadeTime);
end

% Main theme
l = length(mRNAmusicmaintheme);
lhmaintheme = zeros(L*l,1); % array preallocation
for k = 1:l
    lhmaintheme((k-1)*L+1:k*L) = fnote(mRNAmusicmaintheme{k}, t, fs, fadeTime);
end

% Finale
l = length(mRNAmusicfinale);
lhfinale = zeros(L*l,1); % array preallocation
for k = 1:l
    lhfinale((k-1)*L+1:k*L) = fnote(mRNAmusicfinale{k}, t, fs, fadeTime);
end

disp('....Done')

% -------------------------
% Protein sonification
% -------------------------
disp('Protein sonification')
fprintf('....Convert %i amino acids into a notes array\n', length(protein))

% Array preallocation. The array must be 1 cell longer because the last codon is
% a stop codon: so there is not amino acid (insert a rest).
Proteinmusic = cell(length(protein)+1,4);
Proteinmusic(:,1:3) = {'r'};

Proteinmusic(protein=='I',4) = {'B4'};
Proteinmusic(protein=='V',4) = {'A4'};
Proteinmusic(protein=='L',4) = {'F#4'};
Proteinmusic(protein=='M',4) = {'E4'};
Proteinmusic(protein=='F',4) = {'D4'};
Proteinmusic(protein=='W',4) = {'B3'};
Proteinmusic(protein=='Y',4) = {'A3'};
Proteinmusic(protein=='C',4) = {'G3'};
Proteinmusic(protein=='A',4) = {'E3'};
Proteinmusic(protein=='P',4) = {'D3'};
Proteinmusic(protein=='G',4) = {'C3'};
Proteinmusic(protein=='T',4) = {'A2'};
Proteinmusic(protein=='S',4) = {'G2'};
Proteinmusic(protein=='Q',4) = {'F2'};
Proteinmusic(protein=='N',4) = {'D2'};
Proteinmusic(protein=='E',4) = {'C2'};
Proteinmusic(protein=='D',4) = {'A#1'};
Proteinmusic(protein=='H',4) = {'G1'};
Proteinmusic(protein=='K',4) = {'F1'};
Proteinmusic(protein=='R',4) = {'D#1'};

Proteinmusic(end,:) = {'r'};
Proteinmusic = Proteinmusic';
Proteinmusic = Proteinmusic(:);

disp('....Done')
disp('....Convert notes array into frequencies array')

% Convert into tunes
l = length(Proteinmusic);
rhmaintheme = zeros(L*l,1); % array preallocation
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
% Reference note A4 = 440 Hz
a = ['A';' ';'B';'C';' ';'D';' ';'E';'F';' ';'G']; % keyboard
iv = find(a == str(1), 1, 'first');                % find the note index
if isempty(iv)
    y = zeros(size(t));
    return
end

oct = str2double(str(end)) * 12;   % octave offset
alt = (length(str) == 3);          % sharp (#) present
key = iv + alt + oct;              % compute the key number (relative indexing)
f = 440 * (2^(1/12))^(key - 49);   % compute frequency

% Create pitch vector
y = sin((2*pi*f) .* t);

% Fading the vector to avoid clicks
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
%PLAYWITHPROGRESS Play audio with a simple progress bar.
% Uses audioplayer to track CurrentSample.

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
