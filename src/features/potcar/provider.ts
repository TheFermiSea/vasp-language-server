import { Diagnostic } from 'vscode-languageserver/node';
import { TextDocument } from 'vscode-languageserver-textdocument';
import { BaseFeatureProvider } from '../../core/feature-provider';
import { validatePotcar } from './linting';

export class PotcarFeatureProvider extends BaseFeatureProvider {
    async validate(document: TextDocument): Promise<Diagnostic[]> {
        return validatePotcar(document);
    }
}
