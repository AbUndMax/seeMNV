package com.github.abundmax.seemnv.core;

import java.util.BitSet;

public record MNVParticipation(MNVCandidate candidate, BitSet[] participationMasks) {}