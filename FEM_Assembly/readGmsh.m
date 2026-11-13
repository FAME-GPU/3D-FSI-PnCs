function meshData = readGmsh(filename, cellSize, solidTag, fluidTag)
if nargin < 2 || isempty(cellSize)
    cellSize = 1.0;
end
if nargin < 3 || isempty(solidTag)
    solidTag = 2;
end
if nargin < 4 || isempty(fluidTag)
    fluidTag = 1;
end

if isscalar(cellSize)
    cellSize = [cellSize, cellSize, cellSize];
end
cellSize = cellSize(:);
if numel(cellSize) ~= 3 || any(~isfinite(cellSize)) || any(cellSize <= 0)
    error('cellSize must be a positive scalar or a positive 3-vector [Lx,Ly,Lz].');
end

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open Gmsh file: %s', filename);
end
cleaner = onCleanup(@() fclose(fid)); 

meshData.nodes = [];
meshData.elements = [];
meshData.cellIndices = [];

while ~feof(fid)
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end

    if strcmp(strtrim(tline), '$Nodes')
        numNodes = fscanf(fid, '%d', 1);
        nodeData = fscanf(fid, '%f', [4, numNodes]);
        if size(nodeData,2) ~= numNodes
            error('Failed while reading the $Nodes section.');
        end
        meshData.nodes = nodeData(2:4,:);
        fgetl(fid); 
        fgetl(fid); 

    elseif strcmp(strtrim(tline), '$Elements')
        numElements = fscanf(fid, '%d', 1);
        fgetl(fid); 

        tetElements = zeros(4,0);
        cellIdx = zeros(1,0);
        for i = 1:numElements
            lineText = fgetl(fid);
            lineData = sscanf(lineText, '%f').';
            if isempty(lineData)
                continue;
            end

            elmType = lineData(2);
            if elmType == 4 
                numTags = lineData(3);
                if numTags < 1
                    error('A tetrahedron has no physical-volume tag.');
                end
                physicalTag = lineData(4);
                nodeStart = 4 + numTags;
                nodesList = lineData(nodeStart:end);
                if numel(nodesList) ~= 4
                    error('Only four-node tetrahedra are supported.');
                end
                tetElements(:,end+1) = nodesList(:); 
                cellIdx(end+1) = physicalTag;
            end
        end
        meshData.elements = tetElements;
        meshData.cellIndices = cellIdx;
        fgetl(fid); 
    end
end

if isempty(meshData.nodes) || isempty(meshData.elements)
    error('No nodes or four-node tetrahedra were found. Use MSH 2.2 ASCII.');
end

cell_s_idx = find(meshData.cellIndices == solidTag);
cell_f_idx = find(meshData.cellIndices == fluidTag);
meshData.solidTag = solidTag;
meshData.fluidTag = fluidTag;
meshData.cellSize = cellSize;

if isempty(cell_s_idx)
    error('No tetrahedra found for solidTag = %d.', solidTag);
end
if isempty(cell_f_idx)
    error('No tetrahedra found for fluidTag = %d.', fluidTag);
end

cell_s_elements = meshData.elements(:,cell_s_idx);
cell_f_elements = meshData.elements(:,cell_f_idx);
meshData.nodeidx_s = unique(cell_s_elements(:));
meshData.nodeidx_f = unique(cell_f_elements(:));
meshData.interfaceNodes = intersect(meshData.nodeidx_s, meshData.nodeidx_f);

TR_fluid = triangulation(cell_f_elements.', meshData.nodes.');
TR_solid = triangulation(cell_s_elements.', meshData.nodes.');
bf_fluid = sort(freeBoundary(TR_fluid),2);
bf_solid = sort(freeBoundary(TR_solid),2);
interface_faces = intersect(bf_fluid, bf_solid, 'rows');
meshData.interfaceTriangles = interface_faces.';

A = diag(cellSize);
C = A \ meshData.nodes;
tol = 1e-5;
meshData.Faces.left   = find(abs(C(1,:) + 0.5) < tol);
meshData.Faces.right  = find(abs(C(1,:) - 0.5) < tol);
meshData.Faces.front  = find(abs(C(2,:) + 0.5) < tol);
meshData.Faces.back   = find(abs(C(2,:) - 0.5) < tol);
meshData.Faces.bottom = find(abs(C(3,:) + 0.5) < tol);
meshData.Faces.top    = find(abs(C(3,:) - 0.5) < tol);

faceCounts = [numel(meshData.Faces.left), numel(meshData.Faces.right), ...
              numel(meshData.Faces.front), numel(meshData.Faces.back), ...
              numel(meshData.Faces.bottom), numel(meshData.Faces.top)];
if any(faceCounts == 0)
    error('At least one periodic face has no nodes. Check cellSize and centring.');
end

fprintf(['Gmsh mesh loaded: nodes=%d, tetrahedra=%d, solid tets=%d, ' ...
         'fluid tets=%d, interface nodes=%d, interface triangles=%d.\n'], ...
        size(meshData.nodes,2), size(meshData.elements,2), ...
        numel(cell_s_idx), numel(cell_f_idx), ...
        numel(meshData.interfaceNodes), size(meshData.interfaceTriangles,2));

end
