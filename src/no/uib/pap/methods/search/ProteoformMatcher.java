package no.uib.pap.methods.search;

import no.uib.pap.model.Proteoform;

public abstract class ProteoformMatcher {

	public abstract Boolean matches(Proteoform iP, Proteoform rP, Long margin);

	public boolean matches(Long iC, Long rC, Long margin) {
		if (iC != null) {
			if (iC == -1L)
				iC = null;
		}
		if (rC != null) {
			if (rC == -1L)
				rC = null;
		}
		if (iC != null && rC != null) {
			if (iC != rC) {
				if (Math.abs(iC - rC) > margin) {
					return false;
				}
			}
		}
		return true;
	}
}
