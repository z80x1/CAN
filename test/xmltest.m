% Create a sample XML document.
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-07-13
% Created        R O Zhurakivsky 2006-?-?

    docNode = com.mathworks.xml.XMLUtils.createDocument('Workbook');
    docRootNode = docNode.getDocumentElement;
    docRootNode.setAttribute('xmlns','urn:schemas-microsoft-com:office:spreadsheet');
    for i=1:3
	thisWorksheet = docNode.createElement('Worksheet');
%	thisWorksheet.setAttribute('ss:Name','Sheet1');
	thisTable = docNode.createElement('Table');
%	thisTable.setAttribute('ss:StyleID','ta1');
	thisColumn = docNode.createElement('Column');
%	thisColumn.setAttribute('ss:StyleID','Default');
	thisRow = docNode.createElement('Row');
%	thisRow.setAttribute('ss:Height','12.8409');
	thisCell = docNode.createElement('Cell');
	thisData = docNode.createElement('Data');
%	thisData.setAttribute('ss:Type','Number');
	thisData.appendChild(docNode.createTextNode(sprintf('%i',i)));


%	thisWorksheet.appendChild(docNode.createTextNode(sprintf('%i',i)));
	docRootNode.appendChild(thisWorksheet);
	thisWorksheet.appendChild(thisTable);
	thisTable.appendChild(thisColumn);
	thisColumn.appendChild(thisRow);
	thisRow.appendChild(thisCell);
	thisCell.appendChild(thisData);
    end
%    docNode.appendChild(docNode.createComment('this is a comment'));
 
    % Save the sample XML document.
    xmlFileName = ['111.xml']; %#ok
    xmlwrite(xmlFileName,docNode);
    edit(xmlFileName);
