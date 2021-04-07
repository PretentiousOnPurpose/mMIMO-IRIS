function HComp = ChannelCompensation(HHat)

HComp = (conj(HHat))./(abs(HHat) .^ 2);

end