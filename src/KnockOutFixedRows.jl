function KnockOutFixedRows(FixedEls,Tn,Tds,Tss)

	FixedFlag=FixedEls.==1
	NotFixedFlag=FixedEls.!=1
	if any(FixedFlag)
		#Computing for fixed els
		for j=1:length(FixedFlag)
			if FixedFlag[j]==true
				Tn[j]=0.0;
				Tds[j]=0.0;
				Tss[j]=0.0;
			end
		end
	end

	return Tn,Tds,Tss
end