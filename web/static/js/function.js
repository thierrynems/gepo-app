function checkSequence(idField,seq){
	//seq take type of params seq=prot(proteine sequence), seq=gene (gene sequence)
    sequenceValue= document.getElementById(idField).value;
    //alert(sequenceValue);
    if(seq=="prot")sequenceValide=new RegExp("^[ARNDCQEGHILKMFPSTWYVArndcqeghilkmfpstwyv]+$");
    else sequenceValide=new RegExp("^[ATCGatcg]+$");
    sequenceValue=sequenceValue.replace(/[_\s]/g, ''); //replace underscore blank, tabulation, \n to ''
    if(!sequenceValide.test(sequenceValue)){
    	
    }
    document.getElementById(idField).value=sequenceValue;
}
