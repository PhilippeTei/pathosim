# Functions for defining test seeking
import numpy as np

# TODO: Can get_seekers be sped up with numba? Might need to break it down into smaller functions
def get_seekers(candidates, prob): 
        ''' 
        Takes as input a set of potential test-seekers. Determines if they actually will seek a test, probabilistically. 

        Args: 
            candidates  (set): Set of people who satisfy a test-seeking criterion

        Returns: 
            seekers     (set): Set of people who will seek a test
        '''
        # For some criteria, we don't calculate who seeks a test. Example: sx has prob=None, since we already calculate who seeks a test among cs_sx/ncs_sx, and cs_sx/ncs_sx are a partition of sx
        if not(prob is None): 
            test_seeker_uids = np.fromiter(candidates, dtype=int)[np.random.binomial(1, prob, len(candidates)).nonzero()[0]] 
            return set(test_seeker_uids.tolist())

        return set()

def get_all_seekers(eligible_dict, seekprobs):
        '''
        Look at who satisfies the test criteria, and determine which of these people seek a test

        Args: 
            seekprobs        (dict)          : Dictionary containing test seeking probability associated with each criterion
            eligible_dict    (dict)          : Dictionary containing sets of people who meet each criterion

        Returns: 
            test_seekers  (dict)    : Dictionary of test seekers for each criterion
        '''
        test_seekers = dict()

        for criterion in eligible_dict.keys():
            test_seekers[criterion] = get_seekers(eligible_dict[criterion], seekprobs[criterion])

        return test_seekers
