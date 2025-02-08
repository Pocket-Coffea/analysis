from pocket_coffea.utils.configurator import Configurator

class CustomConfigurator(Configurator):
    def __init__(
        self,
        workflow,
        parameters,
        datasets,
        skim,
        preselections,
        categories,
        weights,
        variations,
        variables,
        weights_classes=None,    
        columns=None,
        workflow_options=None,
        save_skimmed_files=None,
        do_postprocessing=True,
        lepton="Electron"
    ):
        super().__init__(
            workflow,
            parameters,
            datasets,
            skim,
            preselections,
            categories,
            weights,
            variations,
            variables,
            weights_classes=None,    
            columns=None,
            workflow_options=None,
            save_skimmed_files=None,
            do_postprocessing=True,
        )
        self.lepton = lepton
    